# Provenance

The configurations in `upload/` were taken from the runs, on
`earthsim-01:/data/slee91/rsf-stepcheck-20260831`, its counterpart of the same name on
`nemex-gpu`, and `nemex-gpu:/data/slee91/rsf-wavefield-20260831`. The two hosts have separate
filesystems and the calculations are split across them.

`SOURCE_COMMIT.txt` and `SUBMODULE_COMMITS.txt` are the campaign's own records. This
directory is not part of the upload.

Ninety-three adaptive aging-law configurations set `control.rsf_dtheta_max = 0.2`. The two
wavefield restarts come from an earlier campaign and omit it. The nine fixed-time-step
shear-box cases do not evaluate the adaptive limiter, and the elastoplastic application
baseline has no RSF; inactive limiter entries were removed from those ten files.

All 101 RSF configurations set
`control.rsf_slip_rate_projection_option = 1`, the total-deviatoric-strain invariant measure in
the pinned source commit. The calculations used this same rate path while it was numbered 3;
the final two-option interface renumbers it to 1 without changing the rate definition.

`truerho.patch` is the one-line change the Appendix E.4 comparison was actually run
against, before `control.mass_scaling_reference_speed` existed. The two configurations in
`decollement/density_floor/` reach the same two behaviours through that setting.
