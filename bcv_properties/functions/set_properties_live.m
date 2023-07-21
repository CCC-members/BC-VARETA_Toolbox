function properties = set_properties_live(properties,gp,sp,ap,cp)

% Updating general params

properties.general_params.bcv_workspace.BCV_input_dir        = gp.BCV_input_dir;
properties.general_params.bcv_workspace.BCV_work_dir         = gp.BCV_work_dir;
properties.general_params.analysis_level.value               = gp.analysis_level;
properties.general_params.run_by_trial.value             = gp.run_by_trial;
properties.general_params.run_frequency_bin.value        = gp.run_frequency_bin;
properties.general_params.system_response.value          = gp.system_response;

% Updating sensor level params
properties.sensor_params.method.value                        = sp.method;

properties.sensor_params.frequencies                         = sp.frequencies;
properties.sensor_params.freq_resol.value                    = sp.freq_resol;
properties.sensor_params.samp_freq.value                     = sp.samp_freq;
properties.sensor_params.max_freq.value                      = sp.max_freq;
properties.sensor_params.freq_gfiltvar.value                 = sp.freq_gfiltvar;
properties.sensor_params.win_order.value                     = sp.win_order;

% Updating activation level params
properties.activation_params.methods                         = ap.actv_methods;
properties.activation_params.IsCurv.value                    = ap.act_params.IsCurv;
properties.activation_params.IsParcel.value                  = ap.act_params.IsParcel;
properties.activation_params.IsNeigh.value                   = ap.act_params.IsNeigh;
properties.activation_params.IsField.value                   = ap.act_params.IsField;
properties.activation_params.aGiri.value                     = ap.act_params.aGiri;
properties.activation_params.aSulc.value                     = ap.act_params.aSulc;
properties.activation_params.bGiri.value                     = ap.act_params.bGiri;
properties.activation_params.bSulc.value                     = ap.act_params.bSulc;
properties.activation_params.regLaplacian.value              = ap.act_params.regLaplacian;

% Updating connectivity level params
properties.connectivity_params.methods                       = cp.conn_methods;
properties.connectivity_params.IsCurv.value                  = cp.conn_params.IsCurv;
properties.connectivity_params.IsNeigh.value                 = cp.conn_params.IsNeigh;
properties.connectivity_params.IsField.value                 = cp.conn_params.IsField;
properties.connectivity_params.aGiri.value                   = cp.conn_params.aGiri;
properties.connectivity_params.aSulc.value                   = cp.conn_params.aSulc;
properties.connectivity_params.bGiri.value                   = cp.conn_params.bGiri;
properties.connectivity_params.bSulc.value                   = cp.conn_params.bSulc;
properties.connectivity_params.axi.value                     = cp.conn_params.axi;
properties.connectivity_params.maxiter_outer.value           = cp.conn_params.maxiter_outer;
properties.connectivity_params.maxiter_inner.value           = cp.conn_params.maxiter_inner;
properties.connectivity_params.ntry.value                    = cp.conn_params.ntry;
properties.connectivity_params.prew.value                    = cp.conn_params.prew;
properties.connectivity_params.penalty.value                 = cp.conn_params.penalty;
properties.connectivity_params.rth1.value                    = cp.conn_params.rth1;
properties.connectivity_params.rth2.value                    = cp.conn_params.rth2;
properties.connectivity_params.eigreg.value                  = cp.conn_params.eigreg;
properties.connectivity_params.regLaplacian.value            = cp.conn_params.regLaplacian;

end

