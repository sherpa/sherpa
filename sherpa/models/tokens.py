#
#  Copyright (C) 2021, 2023
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""Helper routines to reduce the number of brackets in a model expression.

The model expression is a tree, but it's easier (to me) to deal with
a list, so we create a stream of tokens:

  - term (i.e. a single model instance or a numeric term)
  - operator
  - left bracket
  - right bracket.
"""

import numpy as np

__all__ = ('simplify', )


class Token:
    """Represent a model expression as tokens."""

    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

    # __repr__ is for debugging only
    def __repr__(self):
        return f'{self.__class__.__name__}("{self.name}")'


class TokenTerm(Token):
    """A term"""

    pass


class TokenOp(Token):
    """An operator"""

    def __init__(self, name, precedence):
        self.precedence = precedence
        super().__init__(name)

    def __repr__(self):
        return f'{self.__class__.__name__}("{self.name}", {self.precedence})'


class TokenUnOp(TokenOp):
    """A unary operator"""

    pass


class TokenBinOp(TokenOp):
    """A binary operator"""

    pass


class TokenLeft(Token):
    """A left bracket

    We are sent the operator precedence to the left of the bracket
    (may be None) and the operator precedence for the
    exclosed term.

    """

    def __init__(self, ctr, precedence, left=None):
        self.counter = ctr
        self.precedence = precedence
        self.left = left
        super().__init__('(')

    def __repr__(self):
        return f'TokenLeft({self.counter}, {self.precedence}, {self.left})'


class TokenRight(Token):
    """A right bracket

    We are sent the operator precedence for the exclosed term and
    operator precedence to the right of the bracket. It is expected
    that the right operator precedence is not set when the object
    is created.

    """

    def __init__(self, ctr, precedence, right=None):
        self.counter = ctr
        self.precedence = precedence
        self.right = right
        super().__init__(')')

    def __repr__(self):
        return f'TokenRight({self.counter}, {self.precedence}, {self.right})'


class Store:
    """Supply an increasing integer value."""

    def __init__(self):
        self.store = 0

    def get_counter(self):
        """Return the counter and increase the state"""

        ctr = self.store
        self.store += 1
        return ctr


def bracket(store, tokens, precedence, left=None):
    """Add brackets around an expression.

    Parameters
    ----------
    store : Store instance
        Used to create matching brackets.
    tokens : list of Token
        The expression to wrap.
    precedence : int
        The precedence of the expression
    left : int or None, optional
        The precedence to the left of this bracket.
    force : bool, optional
        If True then the bracket is always added.

    Notes
    -----
    This makes a small simplification, as we don't add brackets around
    a term.

    """

    # Is this needed or does it just simplify some things but at the
    # expense of complicating others?
    if len(tokens) == 1 and isinstance(tokens[0], TokenTerm):
        return tokens

    ctr = store.get_counter()
    return [TokenLeft(ctr, precedence, left)] + tokens + [TokenRight(ctr, precedence)]


# This MUST match sherpa.models.model.op_to_precedence(numpy.subtract)
# but it is not written as such to avoid import loops.
#
ZERO_PRECEDENCE = 0


def left_token(store, model, leftprec=None):
    """Tokenize the expression but do not process the right brackets.

    We also handle repeated subtraction terms, where we want the
    precedence to decrease (so a - (b - c) has the 'b - c' term with a
    lower precedence).

    Parameters
    ----------
    store : Store instance
        Used to create matching brackets.
    model : Model instance
    leftprec : int or None, optional
        The precedence to the left of the expression.

    Returns
    -------
    tokens : list of Token
        The converted stream.

    Notes
    -----
    A pass needs to be made to the output to add in the correct "right
    precedence" for TokenRight elements. This is done by right_token.

    """

    # Rather than use instance checks, try to use a more-Pythonic
    # approach. UnaryOp models have a single parts element, BinaryOp
    # have two, and others have no parts element. Alternatively we
    # could have looked for .arg (unary) or .lhs and .rhs (binary).
    # However, there are other composite models with 1 value that do
    # not follow UnaryUpModel, so just access the parts array.
    #
    try:
        parts = model.parts
    except AttributeError:
        return [TokenTerm(model._orig_name)]

    # UnaryOpModel and BinaryOpModel require processing their contents,
    # but other Composite models may not require any analysis (they
    # send in the constituent .name attribute to some model-specific
    # string). We use the presence of the op attribute as a
    # way to say "this is a *OpModel" case.
    #
    if not hasattr(model, 'op'):
        return [TokenTerm(model._orig_name)]

    if len(parts) == 1:
        try:
            oprec = model.precedence
        except AttributeError:
            oprec = None

        out = [TokenUnOp(model.opstr, oprec)]

        # Note: we remove the left-precedence here as it should be
        # "separate" from the existing left-precedence.
        #
        tokens = left_token(store, model.parts[0], None)

        # For functions we require a bracket, so we can just treat this
        # as a term. For operators we would like to remove excess brackets,
        # but the rules for the binary operator don't work here. A
        # simple solution is to only remove brackets when a single term
        # is being operated on.
        #
        # I think that we could just check the length and we do not need
        # the isinstance check, but leave in.
        if model.op in [np.negative, np.positive] and len(tokens) == 1 \
           and isinstance(tokens[0], TokenTerm):
            rhs = tokens
        else:
            rhs = [TokenTerm('(')] + tokens + [TokenTerm(')')]

        return out + rhs

    if len(parts) == 2:
        lhs = model.parts[0]
        rhs = model.parts[1]
        try:
            lprec = lhs.get_precedence()
        except AttributeError:
            lprec = None

        try:
            rprec = rhs.get_precedence()
        except AttributeError:
            rprec = None

        # By construction we should have a precedence, but
        # add a fall-through just in case.
        #
        try:
            oprec = model.precedence
        except AttributeError:
            oprec = None

        # Do we have to deal with subtracting a term which itself
        # contains a subtraction? We use a <= check in case we are in
        # a chain of changes. I am not convinced this is handled
        # correctly.
        #
        # If we have a subtraction of a term which contains a
        # subtraction, we need to process those terms to adjust the
        # precedences. We do this to the token list as we do not want
        # to change the BinaryOpModel in case it is used in multiple
        # expressions.
        #
        ltokens = left_token(store, lhs, leftprec)
        lterms = bracket(store, ltokens, lprec, left=leftprec)

        rtokens = left_token(store, rhs, oprec)

        if oprec is not None and rprec is not None and oprec == ZERO_PRECEDENCE and rprec == ZERO_PRECEDENCE:
            ctr = -1
            for token in rtokens:
                if not hasattr(token, 'precedence'):
                    continue

                if token.name == '-':
                    # Not sure about this
                    if token.precedence <= ctr:
                        ctr = token.precedence - 1

                    token.precedence = ctr
                    ctr -= 1

            rprec = ctr  # Not sure about this

        rterms = bracket(store, rtokens, rprec, left=oprec)

        out = [TokenBinOp(model.opstr, oprec)]
        return lterms + out + rterms

    # We should never get here, but if we do, do not error out, just
    # do nothing.
    #
    return [TokenTerm(model._orig_name)]


def right_token(tokens):
    """Process a list of tokens to add the right-bracket precedences.

    Parameters
    ----------
    tokens : list of Tokens

    Returns
    -------
    tokens : list of Tokens

    """

    # - reverse the list
    # - precedence = None
    # - loop through each element
    #   - if TokenOp, set precedence to this value and copy over
    #   - if TokenRight, update the right precedence
    #   - otherwise copy over
    # - reverse the list
    #
    out = []
    precedence = None
    for tok in tokens[::-1]:

        # Do we want unary-op precedences to pass through? Probably.
        if isinstance(tok, TokenOp):
            precedence = tok.precedence

        if not isinstance(tok, TokenRight):
            out.append(tok)
            continue

        out.append(TokenRight(tok.counter, tok.precedence, right=precedence))

    return out[::-1]


def tokenize(model):
    """Convert a model expression into tokens.

    This converts the tree structure into a list.

    Parameters
    ----------
    model : Model instance
        The model expression

    Returns
    -------
    tokens : list of Token
        The expression converted to tokens

    """

    store = Store()
    return right_token(left_token(store, model))


def simplify_brackets(tokens):
    """Remove excess brackets.

    Parameters
    ----------
    tokens : list of Token
        The output of tokenize.

    Returns
    -------
    expr : string
        The model expression
    """

    # Need to extract the left and right precedences for each
    # bracket pair.
    #
    brackets = {}
    for tok in tokens:
        if isinstance(tok, TokenLeft):
            assert tok.counter not in brackets
            brackets[tok.counter] = {'left': tok.left}
            continue

        if isinstance(tok, TokenRight):
            brackets[tok.counter]['right'] = tok.right
            continue

    # Add each term to the output stream unless it's a bracket,
    # in which case compare the precedence stored within the
    # expression to the left/right operators.
    #
    out = ""
    for tok in tokens:
        if isinstance(tok, TokenTerm):
            out += tok.name
            continue

        if isinstance(tok, TokenUnOp):
            out += f'{tok.name}'
            continue

        if isinstance(tok, TokenBinOp):
            out += f' {tok.name} '
            continue

        assert isinstance(tok, (TokenLeft, TokenRight))

        # I have removed the simplification in bracket which means we
        # now have terms without a precedence.
        #
        if tok.precedence is None:
            # drop the bracket
            continue

        b = brackets[tok.counter]
        lprec = b['left']
        rprec = b['right']

        # We want to check if the precedence within the bracket is
        # less than either the left or right brackets. This is
        # slightly complicated by the presence of None, as we want to
        # ignore these values.
        #
        lnone = lprec is None
        rnone = rprec is None
        if lnone and rnone:
            continue

        if lnone:
            external = rprec
        elif rnone:
            external = lprec
        else:
            external = min(lprec, rprec)

        if tok.precedence < external:
            out += tok.name

    return out


def simplify(expr):
    """Remove excess brackets.

    Parameters
    ----------
    expr : Model instance
        The model.

    Returns
    -------
    expr : string
        The model expression

    Notes
    -----

    A complex expression will call simplify multiple times since each
    sub-element (e.g. BinaryOpModel) will access the model's name
    attribute.  Some composite models, such as RSPModel, create a name
    using a format like

        apply_rmf(apply_arf(xxx.name))

    where xxx is the sub-component, and others like BinaryOpModel use

        (lhs.name op rhs.name)

    """

    # This is not a perfect check, but it should catch obvious problems.
    #
    if isinstance(expr, str):
        raise ValueError(f"simplify sent a string: {expr}")

    toks = tokenize(expr)
    return simplify_brackets(toks)
