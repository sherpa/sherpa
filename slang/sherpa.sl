%
%_SLANG_INSERT_SAO_COPYRIGHT_HERE_(2007)_
%_SLANG_INSERT_GPL_LICENSE_HERE_
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% THIS IS A GENERATED FILE -- DO NOT EDIT!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


require("pysl");

py_import("pysl");  % Allow pickling of S-Lang objects
py_import("sherpa.astro.ui");

static variable _session = py.sherpa.astro.ui._session;

#ifnexists SherpaErr
new_exception("SherpaErr", RunTimeError, "Sherpa Err");
#endif
py2sl_errormap[py.sherpa.utils.err.SherpaErr] = SherpaErr;

#ifnexists SherpaDataErr
new_exception("SherpaDataErr", SherpaErr, "Sherpa Data Err");
#endif
py2sl_errormap[py.sherpa.utils.err.DataErr] = SherpaDataErr;

#ifnexists SherpaIOErr
new_exception("SherpaIOErr", SherpaErr, "Sherpa IO Err");
#endif
py2sl_errormap[py.sherpa.utils.err.IOErr] = SherpaIOErr;

#ifnexists SherpaNotImplementedErr
new_exception("SherpaNotImplementedErr", SherpaErr, "Sherpa NotImplemented Err");
#endif
py2sl_errormap[py.sherpa.utils.err.NotImplementedErr] = SherpaNotImplementedErr;

#ifnexists SherpaArgumentErr
new_exception("SherpaArgumentErr", SherpaErr, "Sherpa Argument Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ArgumentErr] = SherpaArgumentErr;

#ifnexists SherpaArgumentTypeErr
new_exception("SherpaArgumentTypeErr", SherpaErr, "Sherpa ArgumentType Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ArgumentTypeErr] = SherpaArgumentTypeErr;

#ifnexists SherpaConfidenceErr
new_exception("SherpaConfidenceErr", SherpaErr, "Sherpa Confidence Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ConfidenceErr] = SherpaConfidenceErr;

#ifnexists SherpaDS9Err
new_exception("SherpaDS9Err", SherpaErr, "Sherpa DS9 Err");
#endif
py2sl_errormap[py.sherpa.utils.err.DS9Err] = SherpaDS9Err;

#ifnexists SherpaEstErr
new_exception("SherpaEstErr", SherpaErr, "Sherpa Est Err");
#endif
py2sl_errormap[py.sherpa.utils.err.EstErr] = SherpaEstErr;

#ifnexists SherpaFitErr
new_exception("SherpaFitErr", SherpaErr, "Sherpa Fit Err");
#endif
py2sl_errormap[py.sherpa.utils.err.FitErr] = SherpaFitErr;

#ifnexists SherpaIdentifierErr
new_exception("SherpaIdentifierErr", SherpaErr, "Sherpa Identifier Err");
#endif
py2sl_errormap[py.sherpa.utils.err.IdentifierErr] = SherpaIdentifierErr;

#ifnexists SherpaImportErr
new_exception("SherpaImportErr", SherpaErr, "Sherpa Import Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ImportErr] = SherpaImportErr;

#ifnexists SherpaInstrumentErr
new_exception("SherpaInstrumentErr", SherpaErr, "Sherpa Instrument Err");
#endif
py2sl_errormap[py.sherpa.utils.err.InstrumentErr] = SherpaInstrumentErr;

#ifnexists SherpaModelErr
new_exception("SherpaModelErr", SherpaErr, "Sherpa Model Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ModelErr] = SherpaModelErr;

#ifnexists SherpaPSFErr
new_exception("SherpaPSFErr", SherpaErr, "Sherpa PSF Err");
#endif
py2sl_errormap[py.sherpa.utils.err.PSFErr] = SherpaPSFErr;

#ifnexists SherpaParameterErr
new_exception("SherpaParameterErr", SherpaErr, "Sherpa Parameter Err");
#endif
py2sl_errormap[py.sherpa.utils.err.ParameterErr] = SherpaParameterErr;

#ifnexists SherpaPlotErr
new_exception("SherpaPlotErr", SherpaErr, "Sherpa Plot Err");
#endif
py2sl_errormap[py.sherpa.utils.err.PlotErr] = SherpaPlotErr;

#ifnexists SherpaSessionErr
new_exception("SherpaSessionErr", SherpaErr, "Sherpa Session Err");
#endif
py2sl_errormap[py.sherpa.utils.err.SessionErr] = SherpaSessionErr;

#ifnexists SherpaStatErr
new_exception("SherpaStatErr", SherpaErr, "Sherpa Stat Err");
#endif
py2sl_errormap[py.sherpa.utils.err.StatErr] = SherpaStatErr;


#if (_slang_version >= 20101)
%Return all the qualifiers in a struct
private define _sherpa_get_qualifiers (kwstruct);
private define _sherpa_get_qualifiers (kwstruct)
{
 if( kwstruct == NULL )
    return NULL;
  variable ii, name;
  variable kwnames = get_struct_field_names(kwstruct);
  variable kwargs = Struct_Type[2 * length(kwnames)];
  ii = 0;
  foreach (kwnames)
  {
     name = ();
     kwargs[ii] = struct{value=name};
     ii = ii + 1;
     kwargs[ii] = struct{value=get_struct_field(kwstruct, name)};
     ii = ii + 1;
   }
 return kwargs;
}
#endif


variable Data1D = py.sherpa.astro.ui.Data1D;
        
variable Data1DInt = py.sherpa.astro.ui.Data1DInt;
        
variable Data2D = py.sherpa.astro.ui.Data2D;
        
variable Data2DInt = py.sherpa.astro.ui.Data2DInt;
        
variable DataARF = py.sherpa.astro.ui.DataARF;
        
variable DataIMG = py.sherpa.astro.ui.DataIMG;
        
variable DataPHA = py.sherpa.astro.ui.DataPHA;
        
variable DataRMF = py.sherpa.astro.ui.DataRMF;
        
variable Prior = py.sherpa.astro.ui.Prior;
        
define add_user_pars() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: add_user_pars( modelname, parnames, [ parvals, [ parmins, [ parmaxs, [ parunits, [ parfrozen]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.add_user_pars, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.add_user_pars, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable atten = py.sherpa.astro.ui.atten;
        
variable bbody = py.sherpa.astro.ui.bbody;
        
variable bbodyfreq = py.sherpa.astro.ui.bbodyfreq;
        
variable beta1d = py.sherpa.astro.ui.beta1d;
        
variable beta2d = py.sherpa.astro.ui.beta2d;
        
variable box1d = py.sherpa.astro.ui.box1d;
        
variable box2d = py.sherpa.astro.ui.box2d;
        
variable bpl1d = py.sherpa.astro.ui.bpl1d;
        
define calc_chisqr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_chisqr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_chisqr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_data_sum() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_data_sum, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_data_sum, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_data_sum2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_data_sum2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_data_sum2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_energy_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_energy_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_energy_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_ftest() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_ftest, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_ftest, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_kcorr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: calc_kcorr( z, obslo, obshi, [ restlo, [ resthi, [ id, [ bkg_id]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_kcorr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_kcorr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_mlr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_mlr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_mlr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_model_sum() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_model_sum, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_model_sum, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_model_sum2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_model_sum2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_model_sum2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_photon_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_photon_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_photon_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_source_sum() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_source_sum, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_source_sum, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_source_sum2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_source_sum2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_source_sum2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define calc_stat() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.calc_stat, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.calc_stat, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define clean() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.clean, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.clean, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define conf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.conf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.conf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define confidence() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.confidence, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.confidence, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable const1d = py.sherpa.astro.ui.const1d;
        
variable const2d = py.sherpa.astro.ui.const2d;
        
define contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_fit_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_fit_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_fit_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_kernel() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_kernel, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_kernel, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_ratio() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_ratio, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_ratio, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define contour_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.contour_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.contour_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define copy_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: copy_data( fromid, toid )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.copy_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.copy_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
static variable cos = py.sherpa.astro.ui.cos;
        
define covar() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.covar, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.covar, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define covariance() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.covariance, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.covariance, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define create_model_component() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.create_model_component, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.create_model_component, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define dataspace1d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: dataspace1d( start, stop, [ step, [ numbins, [ id, [ bkg_id, [ dstype]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.dataspace1d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.dataspace1d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define dataspace2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: dataspace2d( dims, [ id, [ dstype]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.dataspace2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.dataspace2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define delete_bkg_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.delete_bkg_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.delete_bkg_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define delete_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.delete_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.delete_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define delete_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.delete_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.delete_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define delete_model_component() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: delete_model_component( name )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.delete_model_component, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.delete_model_component, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define delete_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.delete_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.delete_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable delta1d = py.sherpa.astro.ui.delta1d;
        
variable delta2d = py.sherpa.astro.ui.delta2d;
        
variable dered = py.sherpa.astro.ui.dered;
        
variable devaucouleurs2d = py.sherpa.astro.ui.devaucouleurs2d;
        
variable edge = py.sherpa.astro.ui.edge;
        
define eqwidth() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: eqwidth( src, combo, [ id, [ lo, [ hi, [ bkg_id]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.eqwidth, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.eqwidth, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable erf = py.sherpa.astro.ui.erf;
        
variable erfc = py.sherpa.astro.ui.erfc;
        
static variable exp = py.sherpa.astro.ui.exp;
        
variable exp10 = py.sherpa.astro.ui.exp10;
        
define fake() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.fake, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.fake, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define fake_pha() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: fake_pha( id, arf, rmf, exposure, [ backscal, [ areascal, [ grouping, [ grouped, [ quality, [ bkg]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.fake_pha, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.fake_pha, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define fit_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.fit_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.fit_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define freeze() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.freeze, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.freeze, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define gamma() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.gamma, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.gamma, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable gauss1d = py.sherpa.astro.ui.gauss1d;
        
variable gauss2d = py.sherpa.astro.ui.gauss2d;
        
define get_analysis() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_analysis, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_analysis, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_areascal() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_areascal, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_areascal, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_arf_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_arf_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_arf_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_axes() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_axes, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_axes, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_backscal() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_backscal, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_backscal, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_chisqr_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_chisqr_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_chisqr_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_delchi_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_delchi_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_delchi_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_fit_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_fit_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_fit_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_model_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_model_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_model_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_ratio_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_ratio_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_ratio_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_resid_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_resid_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_resid_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_bkg_source_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_bkg_source_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_bkg_source_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_chisqr_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_chisqr_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_chisqr_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_conf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_conf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_conf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_conf_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_conf_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_conf_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_conf_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_conf_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_conf_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_confidence_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_confidence_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_confidence_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_coord() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_coord, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_coord, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_counts() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_counts, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_counts, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_covar() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_covar, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_covar, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_covar_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_covar_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_covar_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_covar_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_covar_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_covar_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_covariance_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_covariance_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_covariance_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data_contour_prefs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data_contour_prefs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data_contour_prefs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_data_plot_prefs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_data_plot_prefs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_data_plot_prefs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_default_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_default_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_default_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_delchi_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_delchi_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_delchi_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_dep() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_dep, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_dep, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_dims() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_dims, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_dims, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_energy_flux_hist() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_energy_flux_hist, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_energy_flux_hist, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_error() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_error, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_error, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_exposure() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_exposure, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_exposure, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_filter() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_filter, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_filter, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_fit_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_fit_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_fit_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_fit_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_fit_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_fit_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_fit_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_fit_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_fit_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_functions() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_functions, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_functions, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_grouping() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: get_grouping( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_grouping, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_grouping, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_indep() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_indep, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_indep, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_int_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_int_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_int_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_int_unc() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_int_unc, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_int_unc, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_kernel_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_kernel_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_kernel_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_kernel_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_kernel_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_kernel_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_kernel_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_kernel_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_kernel_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_method() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_method, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_method, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_method_name() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_method_name, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_method_name, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_method_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_method_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_method_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_autoassign_func() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_autoassign_func, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_autoassign_func, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_contour_prefs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_contour_prefs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_contour_prefs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_pars() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: get_model_pars( model )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_pars, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_pars, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_plot_prefs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_plot_prefs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_plot_prefs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_model_type() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: get_model_type( model )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_model_type, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_model_type, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_num_par() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_num_par, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_num_par, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_num_par_frozen() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_num_par_frozen, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_num_par_frozen, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_num_par_thawed() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_num_par_thawed, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_num_par_thawed, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_order_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_order_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_order_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_par() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: get_par( par )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_par, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_par, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_photon_flux_hist() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_photon_flux_hist, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_photon_flux_hist, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_pileup_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_pileup_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_pileup_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_proj_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_proj_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_proj_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_proj_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_proj_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_proj_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_projection_results() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_projection_results, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_projection_results, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_psf_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_psf_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_psf_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_psf_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_psf_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_psf_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_psf_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_psf_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_psf_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_quality() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: get_quality( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_quality, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_quality, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_rate() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_rate, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_rate, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_ratio_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_ratio_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_ratio_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_ratio_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_ratio_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_ratio_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_ratio_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_ratio_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_ratio_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_reg_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_reg_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_reg_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_reg_unc() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_reg_unc, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_reg_unc, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_resid_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_resid_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_resid_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_resid_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_resid_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_resid_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_resid_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_resid_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_resid_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_source_contour() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_source_contour, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_source_contour, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_source_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_source_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_source_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_source_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_source_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_source_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_specresp() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_specresp, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_specresp, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_split_plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_split_plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_split_plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_stat() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_stat, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_stat, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_stat_name() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_stat_name, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_stat_name, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_staterror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_staterror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_staterror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_syserror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_syserror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_syserror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_xsabund() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_xsabund, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_xsabund, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_xscosmo() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_xscosmo, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_xscosmo, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define get_xsxsect() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.get_xsxsect, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.get_xsxsect, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_adapt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_adapt( id, [ min, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_adapt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_adapt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_adapt_snr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_adapt_snr( id, [ min, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_adapt_snr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_adapt_snr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_bins() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_bins( id, [ num, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_bins, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_bins, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_counts() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_counts( id, [ num, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_counts, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_counts, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_snr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_snr( id, [ snr, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_snr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_snr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define group_width() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: group_width( id, [ num, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.group_width, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.group_width, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define guess() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.guess, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.guess, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define histogram1d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.histogram1d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.histogram1d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define histogram2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.histogram2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.histogram2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable hubblereynolds = py.sherpa.astro.ui.hubblereynolds;
        
define igam() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.igam, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.igam, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define igamc() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.igamc, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.igamc, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ignore() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ignore, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ignore, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ignore2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ignore2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ignore2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ignore2d_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: ignore2d_id( ids, [ val] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ignore2d_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ignore2d_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ignore_bad() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ignore_bad, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ignore_bad, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ignore_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: ignore_id( ids, [ lo, [ hi, [ **kwargs]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ignore_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ignore_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_close() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_close, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_close, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_deleteframes() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_deleteframes, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_deleteframes, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_getregion() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_getregion, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_getregion, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_kernel() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_kernel, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_kernel, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_open() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_open, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_open, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_ratio() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_ratio, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_ratio, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_setregion() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: image_setregion( reg )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_setregion, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_setregion, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_xpaget() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: image_xpaget( arg )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_xpaget, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_xpaget, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define image_xpaset() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: image_xpaset( arg )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.image_xpaset, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.image_xpaset, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define incbet() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.incbet, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.incbet, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define int_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: int_proj( par, [ id, [ otherids, [ replot, [ fast, [ min, [ max, [ nloop, [ delv, [ fac, [ log, [ overplot]]]]]]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.int_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.int_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define int_unc() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: int_unc( par, [ id, [ otherids, [ replot, [ min, [ max, [ nloop, [ delv, [ fac, [ log, [ overplot]]]]]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.int_unc, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.int_unc, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable jdpileup = py.sherpa.astro.ui.jdpileup;
        
define lgam() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.lgam, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.lgam, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable linebroad = py.sherpa.astro.ui.linebroad;
        
define link() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: link( par, val )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.link, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.link, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_bkg_ids() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_bkg_ids, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_bkg_ids, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_data_ids() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_data_ids, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_data_ids, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_functions() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_functions, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_functions, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_methods() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_methods, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_methods, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_model_components() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_model_components, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_model_components, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_model_ids() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_model_ids, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_model_ids, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_models() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_models, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_models, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_response_ids() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_response_ids, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_response_ids, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define list_stats() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.list_stats, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.list_stats, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_arf( id, [ arg, [ resp_id, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_arrays() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_arrays( id, [ *args] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_arrays, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_arrays, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_ascii() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_ascii( id, [ filename, [ ncols, [ colkeys, [ dstype, [ sep, [ comment]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_ascii, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_ascii, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_bkg( id, [ arg, [ use_errors, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_bkg_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_bkg_arf( id, [ arg] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_bkg_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_bkg_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_bkg_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_bkg_rmf( id, [ arg] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_bkg_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_bkg_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_conv() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_conv( modelname, filename_or_model, [ *args, [ **kwargs]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_conv, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_conv, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_data( id, [ filename, [ *args, [ **kwargs]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_filter() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_filter( id, [ filename, [ bkg_id, [ *args, [ **kwargs]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_filter, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_filter, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_grouping() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_grouping( id, [ filename, [ bkg_id, [ *args, [ **kwargs]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_grouping, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_grouping, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_image( id, [ arg, [ coord]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_multi_arfs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_multi_arfs( id, filenames, [ resp_ids] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_multi_arfs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_multi_arfs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_multi_rmfs() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_multi_rmfs( id, filenames, [ resp_ids] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_multi_rmfs, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_multi_rmfs, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_pha() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_pha( id, [ arg, [ use_errors]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_pha, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_pha, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_psf( modelname, filename_or_model, [ *args, [ **kwargs]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_quality() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_quality( id, [ filename, [ bkg_id, [ *args, [ **kwargs]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_quality, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_quality, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_rmf( id, [ arg, [ resp_id, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_staterror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_staterror( id, [ filename, [ bkg_id, [ *args, [ **kwargs]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_staterror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_staterror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_syserror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_syserror( id, [ filename, [ bkg_id, [ *args, [ **kwargs]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_syserror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_syserror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_table() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_table( id, [ filename, [ ncols, [ colkeys, [ dstype]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_table, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_table, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_table_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_table_model( modelname, filename, [ *args, [ **kwargs]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_table_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_table_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_user_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_user_model( func, modelname, [ filename, [ *args, [ **kwargs]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_user_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_user_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define load_user_stat() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: load_user_stat( statname, calc_stat_func, [ calc_err_func, [ priors]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.load_user_stat, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.load_user_stat, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
static variable log = py.sherpa.astro.ui.log;
        
static variable log10 = py.sherpa.astro.ui.log10;
        
variable lorentz1d = py.sherpa.astro.ui.lorentz1d;
        
variable lorentz2d = py.sherpa.astro.ui.lorentz2d;
        
variable multiresponsesummodel = py.sherpa.astro.ui.multiresponsesummodel;
        
variable normbeta1d = py.sherpa.astro.ui.normbeta1d;
        
variable normgauss1d = py.sherpa.astro.ui.normgauss1d;
        
define notice() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.notice, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.notice, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define notice2d() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.notice2d, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.notice2d, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define notice2d_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: notice2d_id( ids, [ val] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.notice2d_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.notice2d_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define notice_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: notice_id( ids, [ lo, [ hi, [ **kwargs]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.notice_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.notice_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define pack_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.pack_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.pack_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define pack_pha() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.pack_pha, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.pack_pha, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define pack_table() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.pack_table, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.pack_table, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define paramprompt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.paramprompt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.paramprompt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_chisqr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_chisqr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_chisqr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_delchi() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_delchi, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_delchi, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_fit_delchi() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_fit_delchi, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_fit_delchi, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_fit_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_fit_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_fit_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_ratio() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_ratio, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_ratio, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_bkg_unconvolved() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_bkg_unconvolved, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_bkg_unconvolved, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_chisqr() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_chisqr, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_chisqr, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_delchi() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_delchi, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_delchi, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_energy_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_energy_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_energy_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_fit_delchi() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_fit_delchi, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_fit_delchi, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_fit_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_fit_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_fit_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_kernel() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_kernel, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_kernel, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_order() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_order, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_order, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_photon_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_photon_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_photon_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_ratio() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_ratio, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_ratio, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define plot_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.plot_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.plot_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable poisson = py.sherpa.astro.ui.poisson;
        
variable polynom1d = py.sherpa.astro.ui.polynom1d;
        
variable polynom2d = py.sherpa.astro.ui.polynom2d;
        
variable powlaw1d = py.sherpa.astro.ui.powlaw1d;
        
define proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define projection() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.projection, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.projection, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define rebin() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.rebin, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.rebin, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define reg_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: reg_proj( par0, par1, [ id, [ otherids, [ replot, [ fast, [ min, [ max, [ nloop, [ delv, [ fac, [ log, [ sigma, [ levels, [ overplot]]]]]]]]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.reg_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.reg_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define reg_unc() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: reg_unc( par0, par1, [ id, [ otherids, [ replot, [ min, [ max, [ nloop, [ delv, [ fac, [ log, [ sigma, [ levels, [ overplot]]]]]]]]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.reg_unc, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.reg_unc, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define reset() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: reset( model )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.reset, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.reset, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define restore() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.restore, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.restore, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define sample_energy_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.sample_energy_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.sample_energy_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define sample_photon_flux() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.sample_photon_flux, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.sample_photon_flux, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_all() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_all, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_all, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_arrays() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_arrays( filename, args, [ fields, [ ascii, [ clobber]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_arrays, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_arrays, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_data( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_delchi() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_delchi( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_delchi, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_delchi, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_error() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_error( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_error, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_error, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_filter() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_filter( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_filter, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_filter, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_grouping() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_grouping( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_grouping, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_grouping, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_image( id, [ filename, [ ascii, [ clobber]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_model( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_pha() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_pha( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_pha, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_pha, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_quality() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_quality( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_quality, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_quality, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_resid() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_resid( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_resid, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_resid, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_source( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_staterror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_staterror( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_staterror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_staterror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_syserror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_syserror( id, [ filename, [ bkg_id, [ ascii, [ clobber]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_syserror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_syserror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define save_table() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: save_table( id, [ filename, [ ascii, [ clobber]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.save_table, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.save_table, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable schechter = py.sherpa.astro.ui.schechter;
        
define set_analysis() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_analysis( id, [ quantity, [ type, [ factor]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_analysis, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_analysis, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_areascal() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_areascal( id, [ area, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_areascal, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_areascal, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_arf( id, [ arf, [ resp_id, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_backscal() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_backscal( id, [ backscale, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_backscal, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_backscal, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_bkg( id, [ bkg, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_bkg_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_bkg_model( id, [ model, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_bkg_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_bkg_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_conf_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_conf_opt( name, val )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_conf_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_conf_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_coord() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_coord( id, [ coord] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_coord, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_coord, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_counts() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_dep( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_counts, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_counts, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_covar_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_covar_opt( name, val )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_covar_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_covar_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_data( id, [ data] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_default_id() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_default_id( id )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_default_id, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_default_id, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_dep() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_dep( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_dep, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_dep, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_exposure() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_exposure( id, [ exptime, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_exposure, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_exposure, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_filter() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_filter( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_filter, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_filter, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_grouping() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_grouping( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_grouping, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_grouping, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_method() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_method( meth )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_method, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_method, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_method_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_method_opt( optname, val )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_method_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_method_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_model( id, [ model] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_model_autoassign_func() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_model_autoassign_func, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_model_autoassign_func, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_par() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_par( par, [ val, [ min, [ max, [ frozen]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_par, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_par, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_pileup_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_pileup_model( id, [ model] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_pileup_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_pileup_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_proj_opt() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_proj_opt( name, val )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_proj_opt, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_proj_opt, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_psf( id, [ psf] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_quality() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_quality( id, [ val, [ bkg_id]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_quality, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_quality, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_rmf( id, [ rmf, [ resp_id, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_model( id, [ model] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_stat() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_stat( stat )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_stat, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_stat, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_staterror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_staterror( id, [ val, [ fractional, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_staterror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_staterror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_syserror() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: set_syserror( id, [ val, [ fractional, [ bkg_id]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_syserror, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_syserror, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_xsabund() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_xsabund, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_xsabund, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_xscosmo() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_xscosmo, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_xscosmo, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define set_xsxsect() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.set_xsxsect, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.set_xsxsect, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_all() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_all, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_all, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_bkg_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_bkg_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_bkg_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_bkg_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_bkg_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_bkg_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_conf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_conf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_conf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_covar() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_covar, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_covar, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_filter() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_filter, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_filter, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_fit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_fit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_fit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_kernel() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_kernel, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_kernel, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_method() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_method, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_method, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_model() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_model, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_model, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_proj() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_proj, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_proj, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_psf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_psf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_psf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_source() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_source, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_source, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define show_stat() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.show_stat, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.show_stat, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define simulfit() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.simulfit, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.simulfit, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
static variable sin = py.sherpa.astro.ui.sin;
        
static variable sqrt = py.sherpa.astro.ui.sqrt;
        
variable stephi1d = py.sherpa.astro.ui.stephi1d;
        
variable steplo1d = py.sherpa.astro.ui.steplo1d;
        
define subtract() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.subtract, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.subtract, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable tablemodel = py.sherpa.astro.ui.tablemodel;
        
static variable tan = py.sherpa.astro.ui.tan;
        
define thaw() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.thaw, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.thaw, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define ungroup() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.ungroup, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.ungroup, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unlink() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unlink( par )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unlink, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unlink, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_arf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_arf( arg )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_arf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_arf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_arrays() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_arrays, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_arrays, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_ascii() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_ascii( filename, [ ncols, [ colkeys, [ dstype, [ sep, [ comment]]]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_ascii, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_ascii, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_bkg() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_pha( arg, [ use_errors] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_bkg, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_bkg, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_data() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_data( filename, [ *args, [ **kwargs]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_data, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_data, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_image() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_image( arg, [ coord] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_image, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_image, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_pha() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_pha( arg, [ use_errors] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_pha, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_pha, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_rmf() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_rmf( arg )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_rmf, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_rmf, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unpack_table() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    if( _NARGS == 0 ) {
       () = fprintf(stderr, "%s\n", "Usage: unpack_table( filename, [ ncols, [ colkeys, [ dstype]]] )");
       return;
    }
    
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unpack_table, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unpack_table, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
define unsubtract() {
#if (_slang_version >= 20101)
    variable kwargs = _sherpa_get_qualifiers(__qualifiers());
#else
    variable kwargs = NULL;
#endif
    variable args = __pop_args(_NARGS);
        
    % Interface from S-lang keyword args to PySL keyword args
    if( kwargs == NULL ) {
        py_call(py.sherpa.astro.ui.unsubtract, __push_args(args));
        return;
    }
    else {
        py_call(py.sherpa.astro.ui.unsubtract, __push_args(args),
                PY_KW,  __push_args(kwargs));
    }
}
        
variable usermodel = py.sherpa.astro.ui.usermodel;
        
variable xsabsori = py.sherpa.astro.ui.xsabsori;
        
variable xsacisabs = py.sherpa.astro.ui.xsacisabs;
        
variable xsapec = py.sherpa.astro.ui.xsapec;
        
variable xsbapec = py.sherpa.astro.ui.xsbapec;
        
variable xsbbody = py.sherpa.astro.ui.xsbbody;
        
variable xsbbodyrad = py.sherpa.astro.ui.xsbbodyrad;
        
variable xsbexrav = py.sherpa.astro.ui.xsbexrav;
        
variable xsbexriv = py.sherpa.astro.ui.xsbexriv;
        
variable xsbkn2pow = py.sherpa.astro.ui.xsbkn2pow;
        
variable xsbknpower = py.sherpa.astro.ui.xsbknpower;
        
variable xsbmc = py.sherpa.astro.ui.xsbmc;
        
variable xsbremss = py.sherpa.astro.ui.xsbremss;
        
variable xsbvapec = py.sherpa.astro.ui.xsbvapec;
        
variable xsc6mekl = py.sherpa.astro.ui.xsc6mekl;
        
variable xsc6pmekl = py.sherpa.astro.ui.xsc6pmekl;
        
variable xsc6pvmkl = py.sherpa.astro.ui.xsc6pvmkl;
        
variable xsc6vmekl = py.sherpa.astro.ui.xsc6vmekl;
        
variable xscabs = py.sherpa.astro.ui.xscabs;
        
variable xscemekl = py.sherpa.astro.ui.xscemekl;
        
variable xscevmkl = py.sherpa.astro.ui.xscevmkl;
        
variable xscflow = py.sherpa.astro.ui.xscflow;
        
variable xscompbb = py.sherpa.astro.ui.xscompbb;
        
variable xscompls = py.sherpa.astro.ui.xscompls;
        
variable xscompps = py.sherpa.astro.ui.xscompps;
        
variable xscompst = py.sherpa.astro.ui.xscompst;
        
variable xscomptt = py.sherpa.astro.ui.xscomptt;
        
variable xsconstant = py.sherpa.astro.ui.xsconstant;
        
variable xscutoffpl = py.sherpa.astro.ui.xscutoffpl;
        
variable xscyclabs = py.sherpa.astro.ui.xscyclabs;
        
variable xsdisk = py.sherpa.astro.ui.xsdisk;
        
variable xsdiskbb = py.sherpa.astro.ui.xsdiskbb;
        
variable xsdiskir = py.sherpa.astro.ui.xsdiskir;
        
variable xsdiskline = py.sherpa.astro.ui.xsdiskline;
        
variable xsdiskm = py.sherpa.astro.ui.xsdiskm;
        
variable xsdisko = py.sherpa.astro.ui.xsdisko;
        
variable xsdiskpbb = py.sherpa.astro.ui.xsdiskpbb;
        
variable xsdiskpn = py.sherpa.astro.ui.xsdiskpn;
        
variable xsdust = py.sherpa.astro.ui.xsdust;
        
variable xsedge = py.sherpa.astro.ui.xsedge;
        
variable xsequil = py.sherpa.astro.ui.xsequil;
        
variable xsexpabs = py.sherpa.astro.ui.xsexpabs;
        
variable xsexpdec = py.sherpa.astro.ui.xsexpdec;
        
variable xsexpfac = py.sherpa.astro.ui.xsexpfac;
        
variable xsezdiskbb = py.sherpa.astro.ui.xsezdiskbb;
        
variable xsgabs = py.sherpa.astro.ui.xsgabs;
        
variable xsgaussian = py.sherpa.astro.ui.xsgaussian;
        
variable xsgnei = py.sherpa.astro.ui.xsgnei;
        
variable xsgrad = py.sherpa.astro.ui.xsgrad;
        
variable xsgrbm = py.sherpa.astro.ui.xsgrbm;
        
variable xshighecut = py.sherpa.astro.ui.xshighecut;
        
variable xshrefl = py.sherpa.astro.ui.xshrefl;
        
variable xskerrbb = py.sherpa.astro.ui.xskerrbb;
        
variable xskerrd = py.sherpa.astro.ui.xskerrd;
        
variable xskerrdisk = py.sherpa.astro.ui.xskerrdisk;
        
variable xslaor = py.sherpa.astro.ui.xslaor;
        
variable xslaor2 = py.sherpa.astro.ui.xslaor2;
        
variable xslorentz = py.sherpa.astro.ui.xslorentz;
        
variable xsmeka = py.sherpa.astro.ui.xsmeka;
        
variable xsmekal = py.sherpa.astro.ui.xsmekal;
        
variable xsmkcflow = py.sherpa.astro.ui.xsmkcflow;
        
variable xsnei = py.sherpa.astro.ui.xsnei;
        
variable xsnotch = py.sherpa.astro.ui.xsnotch;
        
variable xsnpshock = py.sherpa.astro.ui.xsnpshock;
        
variable xsnsa = py.sherpa.astro.ui.xsnsa;
        
variable xsnsagrav = py.sherpa.astro.ui.xsnsagrav;
        
variable xsnsatmos = py.sherpa.astro.ui.xsnsatmos;
        
variable xsnsmax = py.sherpa.astro.ui.xsnsmax;
        
variable xsnteea = py.sherpa.astro.ui.xsnteea;
        
variable xsnthcomp = py.sherpa.astro.ui.xsnthcomp;
        
variable xspcfabs = py.sherpa.astro.ui.xspcfabs;
        
variable xspegpwrlw = py.sherpa.astro.ui.xspegpwrlw;
        
variable xspexrav = py.sherpa.astro.ui.xspexrav;
        
variable xspexriv = py.sherpa.astro.ui.xspexriv;
        
variable xsphabs = py.sherpa.astro.ui.xsphabs;
        
variable xsplabs = py.sherpa.astro.ui.xsplabs;
        
variable xsplcabs = py.sherpa.astro.ui.xsplcabs;
        
variable xsposm = py.sherpa.astro.ui.xsposm;
        
variable xspowerlaw = py.sherpa.astro.ui.xspowerlaw;
        
variable xspshock = py.sherpa.astro.ui.xspshock;
        
variable xspwab = py.sherpa.astro.ui.xspwab;
        
variable xsraymond = py.sherpa.astro.ui.xsraymond;
        
variable xsredden = py.sherpa.astro.ui.xsredden;
        
variable xsredge = py.sherpa.astro.ui.xsredge;
        
variable xsrefsch = py.sherpa.astro.ui.xsrefsch;
        
variable xssedov = py.sherpa.astro.ui.xssedov;
        
variable xssmedge = py.sherpa.astro.ui.xssmedge;
        
variable xsspexpcut = py.sherpa.astro.ui.xsspexpcut;
        
variable xsspline = py.sherpa.astro.ui.xsspline;
        
variable xssrcut = py.sherpa.astro.ui.xssrcut;
        
variable xssresc = py.sherpa.astro.ui.xssresc;
        
variable xssss_ice = py.sherpa.astro.ui.xssss_ice;
        
variable xsstep = py.sherpa.astro.ui.xsstep;
        
variable xsswind1 = py.sherpa.astro.ui.xsswind1;
        
variable xstbabs = py.sherpa.astro.ui.xstbabs;
        
variable xstbgrain = py.sherpa.astro.ui.xstbgrain;
        
variable xstbvarabs = py.sherpa.astro.ui.xstbvarabs;
        
variable xsuvred = py.sherpa.astro.ui.xsuvred;
        
variable xsvapec = py.sherpa.astro.ui.xsvapec;
        
variable xsvarabs = py.sherpa.astro.ui.xsvarabs;
        
variable xsvbremss = py.sherpa.astro.ui.xsvbremss;
        
variable xsvequil = py.sherpa.astro.ui.xsvequil;
        
variable xsvgnei = py.sherpa.astro.ui.xsvgnei;
        
variable xsvmcflow = py.sherpa.astro.ui.xsvmcflow;
        
variable xsvmeka = py.sherpa.astro.ui.xsvmeka;
        
variable xsvmekal = py.sherpa.astro.ui.xsvmekal;
        
variable xsvnei = py.sherpa.astro.ui.xsvnei;
        
variable xsvnpshock = py.sherpa.astro.ui.xsvnpshock;
        
variable xsvphabs = py.sherpa.astro.ui.xsvphabs;
        
variable xsvpshock = py.sherpa.astro.ui.xsvpshock;
        
variable xsvraymond = py.sherpa.astro.ui.xsvraymond;
        
variable xsvsedov = py.sherpa.astro.ui.xsvsedov;
        
variable xswabs = py.sherpa.astro.ui.xswabs;
        
variable xswndabs = py.sherpa.astro.ui.xswndabs;
        
variable xsxion = py.sherpa.astro.ui.xsxion;
        
variable xszbbody = py.sherpa.astro.ui.xszbbody;
        
variable xszbremss = py.sherpa.astro.ui.xszbremss;
        
variable xszdust = py.sherpa.astro.ui.xszdust;
        
variable xszedge = py.sherpa.astro.ui.xszedge;
        
variable xszgauss = py.sherpa.astro.ui.xszgauss;
        
variable xszhighect = py.sherpa.astro.ui.xszhighect;
        
variable xszpcfabs = py.sherpa.astro.ui.xszpcfabs;
        
variable xszphabs = py.sherpa.astro.ui.xszphabs;
        
variable xszpowerlw = py.sherpa.astro.ui.xszpowerlw;
        
variable xszredden = py.sherpa.astro.ui.xszredden;
        
variable xszsmdust = py.sherpa.astro.ui.xszsmdust;
        
variable xsztbabs = py.sherpa.astro.ui.xsztbabs;
        
variable xszvarabs = py.sherpa.astro.ui.xszvarabs;
        
variable xszvfeabs = py.sherpa.astro.ui.xszvfeabs;
        
variable xszvphabs = py.sherpa.astro.ui.xszvphabs;
        
variable xszwabs = py.sherpa.astro.ui.xszwabs;
        
variable xszwndabs = py.sherpa.astro.ui.xszwndabs;
        
variable xszxipcf = py.sherpa.astro.ui.xszxipcf;
        
% The autoassign function must be declared public (i.e. global) so that it
% can be pickled by save()
public define _sherpa_assign_model_to_global(name, model) {
    model;
    eval(sprintf("variable %s = ();", name), "Global");
    model.name = sprintf("%s.%s",
                         strlow(py_call(py.type, model).__name__),
                         name);
}

set_model_autoassign_func(&_sherpa_assign_model_to_global);

% create access to help file
$1 = path_concat (path_dirname (__FILE__), "help/sherpa.hlp");
if (NULL != stat_file ($1))
  add_doc_file ($1);
  
variable _sherpa_version = 40201;
variable _sherpa_version_string = "4.2.1";

#ifexists provide
provide("sherpa");
#endif
