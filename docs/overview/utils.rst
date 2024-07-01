***********************
The sherpa.utils module
***********************

.. currentmodule:: sherpa.utils

.. automodule:: sherpa.utils

.. versionchanged:: 4.17.0
   The parameter guess routines, such as `guess_amplitude` and
   `param_apply_limits` should now be taken from the
   :py:mod:`sherpa.utils.guess` module.

.. versionchanged:: 4.16.0
   The `parallel_map` function should now be taken from the
   :py:mod:`sherpa.utils.parallel` module and the numeric types
   (`SherpaFloat`, `SherpaInt`, and `SherpaUInt`) from the
   `sherpa.utils.numeric_types` module.

   .. rubric:: Functions

   .. autosummary::
      :toctree: api

      Knuth_close
      _guess_ampl_scale
      apache_muller
      bisection
      bool_cast
      calc_ftest
      calc_mlr
      calc_total_error
      create_expr
      dataspace1d
      dataspace2d
      demuller
      erf
      export_method
      extract_kernel
      filter_bins
      gamma
      get_error_estimates
      get_fwhm
      get_keyword_defaults
      get_keyword_names
      get_midpoint
      get_num_args
      get_peak
      get_position
      get_valley
      guess_amplitude
      guess_amplitude2d
      guess_amplitude_at_ref
      guess_bounds
      guess_fwhm
      guess_position
      guess_radius
      guess_reference
      histogram1d
      histogram2d
      igam
      igamc
      incbet
      interpolate
      is_binary_file
      lgam
      linear_interp
      multinormal_pdf
      multit_pdf
      nearest_interp
      neville
      neville2d
      new_muller
      normalize
      pad_bounding_box
      parallel_map
      param_apply_limits
      parse_expr
      poisson_noise
      print_fields
      quantile
      rebin
      sao_arange
      sao_fcmp
      set_origin
      sum_intervals
      zeroin

   .. rubric:: Classes

   .. autosummary::
      :toctree: api

      NoNewAttributesAfterInit
