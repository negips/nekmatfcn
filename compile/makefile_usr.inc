SRC_INC := ./src_inc
SRC_USR := ./src_usr
SRC_STATS := ./stats
SRC_ARN := ./src_arn
$(OBJDIR)/rt_filter.o     :$(SRC_USR)/rt_filter.f;                      $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/uparam.o     :$(SRC_USR)/uparam.f;                      $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/IO_tools.o     :$(SRC_USR)/IO_tools.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/chk_div.o     :$(SRC_USR)/chk_div.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/overintegration.o     :$(SRC_USR)/overintegration.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/misc.o     :$(SRC_USR)/misc.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/wing_mvmsh.o     :$(SRC_USR)/wing_mvmsh.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/paul_fsi.o     :$(SRC_USR)/paul_fsi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/chkpoint.o     :$(SRC_USR)/chkpoint.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_2D.o        :$(SRC_STATS)/statistics_2D.f;                         $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_2DIO.o      :$(SRC_STATS)/statistics_2DIO.f;                       $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_2D_usr.o    :$(SRC_STATS)/statistics_2D_usr.f;                     $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_2D_debug.o    :$(SRC_STATS)/statistics_2D_debug.f;                     $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_surf.o    :$(SRC_STATS)/statistics_surf.f;                     $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/statistics_pdf.o    :$(SRC_STATS)/statistics_pdf.f;                     $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/time_series.o          :$(SRC_STATS)/time_series.f;                           $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/time_seriesIO.o        :$(SRC_STATS)/time_seriesIO.f;                         $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/pts_redistribute.o     :$(SRC_STATS)/pts_redistribute.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/uzawa.o     :./$(SRC_USR)/uzawa.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fix_gll.o     :./$(SRC_USR)/fix_gll.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/nl_fsi.o     :./$(SRC_USR)/nl_fsi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/nl_struc.o     :./$(SRC_USR)/nl_struc.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/nl_fsi_pert.o     :./$(SRC_USR)/nl_fsi_pert.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/nl_fsi_arnoldi.o     :./$(SRC_USR)/nl_fsi_arnoldi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/nl_fsi_adjoint.o     :./$(SRC_USR)/nl_fsi_adjoint.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi.o     :$(SRC_USR)/fsi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi_plan3.o     :$(SRC_USR)/fsi_plan3.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi_perturb.o     :$(SRC_USR)/fsi_perturb.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/fsi_arnoldi.o     :$(SRC_USR)/fsi_arnoldi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/err_est.o     :$(SRC_USR)/err_est.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/sfd.o     :$(SRC_USR)/sfd.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/boostconv.o     :$(SRC_USR)/boostconv.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/time_stepper.o     :$(SRC_ARN)/time_stepper.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/arnoldi_arpack.o     :$(SRC_ARN)/arnoldi_arpack.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/conht_tools.o     :$(SRC_ARN)/conht_tools.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/sponge.o     :$(SRC_ARN)/sponge.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/adjnbc.o     :$(SRC_USR)/adjnbc.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/psn_arnoldi.o     :$(SRC_USR)/psn_arnoldi.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/extract_wall_dx.o     :$(SRC_USR)/extract_wall_dx.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/wrapper_lapack.o     :./wrapper_lapack.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@
$(OBJDIR)/matrix_functions.o     :./matrix_functions.f;                      $(F77) -c $(FL2) -I$(SRC_INC) $< -o $@



