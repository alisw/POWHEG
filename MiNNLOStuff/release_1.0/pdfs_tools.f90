!========================================================
!--------------------------------------------------------
! Includes a module providing the parton luminosities
! and their convolutions with splitting and coefficient
! functions @ LO, NLO, NNLO
! -------------------------------------------------------
!========================================================
module pdfs_tools
  use types; use consts_dp
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use rad_tools
  use internal_parameters

  implicit none

  type(dglap_holder), save, public :: dglap_h ! splitting function holder
  type(pdf_table),  save, public :: PDFs  ! parton densities
  type(grid_def), save, public :: grid
  ! here define the building blocks for C1_plus_P
  type(split_mat), save, public :: C1_matrix, P0_matrix
  ! here define the building blocks for C2_plus_P
  type(split_mat), save, public :: C2_matrix, G1_matrix, P1_matrix, P2_matrix
  type(split_mat), save, public :: P0_x_P0_matrix, P0_x_P1_matrix, P1_x_P0_matrix
  type(split_mat), save, public :: C1_x_P0_matrix, C1_x_P1_matrix
  ! now define the various final convolutions
  type(pdf_table),  save, public :: P0_PDFs, P1_PDFs, P2_PDFs
  type(pdf_table),  save, public :: C1_P0_matrix_PDFs, C1_P1_matrix_PDFs
  type(pdf_table),  save, public :: C1_P0_P0_matrix_PDFs
  type(pdf_table),  save, public :: C2_P0_matrix_PDFs
  type(pdf_table),  save, public :: P0_P0_PDFs, P0_P1_PDFs, P1_P0_PDFs
  type(pdf_table),  save, public :: P0_P0_P0_PDFs
  type(pdf_table),  save, public :: C1_P0_PDFs, C1_P1_PDFs
  type(pdf_table),  save, public :: C2_P0_PDFs
  type(pdf_table),  save, public :: G1_P0_PDFs, G1_P1_PDFs
  type(pdf_table),  save, public :: C1_matrix_PDFs, C2_matrix_PDFs
  type(pdf_table),  save, public :: C1_PDFs, C2_PDFs, G1_PDFs ! parton densities convoluted with coefficient functions
  
  private 
  public :: init_pdfs_from_LHAPDF, init_dummy_pdfs, init_pdfs_5FS_evolution, init_pdfs_NNLOPS
  public :: init_C_matxs, reset_dynamical_C_matxs
  public :: get_pdfs_Born, lumi_Born, Dterms, fun_test
  public :: alphas_hoppet, xfxq2_hoppet
  !----------------------------------------------------------------
  real(dp) :: PDF_cutoff
  integer, parameter  :: pdf_fit_order = -5
  real(dp), parameter :: dy = 0.10_dp, ymax = 15.0_dp, dlnlnQ = dy/4.0_dp
  
contains

  ! test function for integration routine
  function fun_test(Ln, x1, x2, msqB, msqV1, msqV2) result(res)
    use types; use consts_dp
    real(dp),                  intent(in) :: Ln, x1, x2
    real(dp), intent(in) :: msqB(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl), msqV1(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl), msqV2(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl)
    real(dp) :: res

    res = -exp(-Ln**2)*2*Ln

    return
  end function fun_test


  !=======================================================================================
  ! set up Hoppet's x grid, splitting functions and PDF tabulation
  ! object (just the memory allocation -- it gets filled below).
  subroutine init_grid_and_dglap(Qmax)
    use coefficient_functions_VH
    real(dp), intent(in), optional :: Qmax
    !-----------------------------------------------------
    type(grid_def) ::  gdarray(4) ! grid

    ! build the PDF grid 
    call InitGridDef(gdarray(4),dy/27.0_dp, 0.2_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(3),dy/9.0_dp,  0.5_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(2),dy/3.0_dp,  2.0_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(1),dy,         ymax  , order=pdf_fit_order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)

    call AllocSplitMat(grid, C1_matrix, nf_lcl)
    call AllocSplitMat(grid, C2_matrix, nf_lcl)
    call AllocSplitMat(grid, G1_matrix, nf_lcl)

    !>> initialise splitting functions
    call qcd_SetNf(nf_lcl)
    call InitDglapHolder(grid, dglap_h,  factscheme = factscheme_MSbar, nloop = 3)

    call InitCoeffMatrix(grid, C1_matrix, C2_matrix, G1_matrix)

    if (present(Qmax)) then
       call AllocPdfTable(grid, PDFs,     PDF_cutoff, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    else
       call AllocPdfTable(grid, PDFs,     PDF_cutoff, 2d4, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    end if
  end subroutine init_grid_and_dglap

  !=======================================================================================
  !! initialise Hoppet's PDF table from LHAPDF and merge it with a hoppet evolution at low Q
  subroutine init_pdfs_NNLOPS(pdf_name, pdf_set, cutoff2_low, cutoff2_high, asmz)
    character(len=*), intent(in)   :: pdf_name
    integer,          intent(in)   :: pdf_set
    real(dp), optional, intent(in) :: asmz, cutoff2_low, cutoff2_high
    !--------------------------------------------------
    !--------------------------------------------------
    type(pdf_table) :: PDFs_lowQ
    real(dp), pointer :: pdf_Q(:,:)
    real(dp) :: alphasMZ, alphasPDF, LHAPDF_cutoff2, LHAPDF_cutoff
    real(dp) :: ref_Q, alphas_Q, threshold_Q
    real(dp), pointer :: pdf_at_ref_Q(:,:)
    integer  :: nf_active, i
    real(dp) :: Qmax
    !--------------------------------------------------
    !>> for tests:
    real(dp) :: xf(-6:6), xtest, Qtest 
    integer  :: itest
    !--------------------------------------------------
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    ! set the PDF cutoff at a low scale
    if (present(cutoff2_low)) then
       PDF_cutoff = sqrt(cutoff2_low)
    else
       PDF_cutoff = 0.6_dp
    end if
    
    !------------------------------------------------------------
    ! set up LHAPDF
    call InitPDFsetByName(trim(pdf_name))
    call InitPDF(pdf_set)

    ! read PDF cutoff from LHAPDF and initialise splitting functions
    if (present(cutoff2_high)) then
       LHAPDF_cutoff = sqrt(cutoff2_high)
    else
       call GetQ2min(pdf_set,LHAPDF_cutoff2)
       LHAPDF_cutoff = sqrt(LHAPDF_cutoff2) * 1.2d0
    end if

    Qmax = max(2e4_dp, cs%rts)
    call init_grid_and_dglap(Qmax)

    ! sort out the coupling: initialise with variable number of active flavours
    if (present(asmz)) then
       alphasMZ = asmz
    else
       alphasMZ = alphasPDF(MZ)
    endif
    call InitRunningCoupling(alfas=alphasMZ, Q=MZ, nloop=3)

    ! fill the PDF grid in the large momentum region:
    call FillPdfTable_LHAPDF(PDFs, evolvePDF)

    !>> compare against LHAPDF
    !!xtest = 0.5_dp
    !!write(*,*) '# test at x = ', xtest
    !!write(*,*) '# LHAPDF:'
    !!do itest = 1, 100
    !!   Qtest = itest/20._dp
    !!   call evolvePDF(xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!write(*,*) ' '
    !!write(*,*) ' '
    !!write(*,*) '# hoppet (LHAPDF evolution):'
    !!do itest = 1, 100
    !!   Qtest = itest/20._dp
    !!   call EvalPdfTable_xQ(PDFs,xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do

    ! now evolve backwards from LHAPDF_cutoff down to PDF_cutoff:    
    ! First check that the number of active flavours is what we expect
    nf_active = 0
    do i = 1, 6
       call GetThreshold(i, threshold_Q)
       if (threshold_Q < LHAPDF_cutoff) nf_active = nf_active + 1
    end do
    if (nf_active .ne. NfAtQ(ash_global, LHAPDF_cutoff)) then
       write(*,*) 'Error pdfs_tools: number of active flavours inconsistent with the PDF thresholds'
       stop
    else
       write(*,*) 'Performing backward evolution from scale Q = ', LHAPDF_cutoff, ', with nf = ', nf_active
    end if
    
    ! extract the PDF at our reference Q value
    ref_Q  = LHAPDF_cutoff
    alphas_Q = alphasPDF(ref_Q)
    call AllocPDF(PDFs%grid, pdf_at_ref_Q)
    call InitPDF_LHAPDF(grid, pdf_at_ref_Q, evolvePDF, ref_Q)

    ! now evolve them and store the result into PDFs_lowQ
    call AllocPdfTable(grid, PDFs_lowQ, PDF_cutoff, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call EvolvePdfTable(PDFs_lowQ, ref_Q, pdf_at_ref_Q, dglap_h, ash_global)

    !!write(*,*) ' '
    !!write(*,*) ' '
    !!write(*,*) '# hoppet (hoppet evolution):'
    !!do itest = 1, 100
    !!   Qtest = itest/20._dp
    !!   call EvalPdfTable_xQ(PDFs_lowQ,xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    
    ! Finally, merge the two PDF grids:
    do i = 0, PDFs%nQ
       ! sanity check
       if (PDFs%Q_vals(i) .ne. PDFs_lowQ%Q_vals(i)) then
          write(*,*) 'Error pdfs_tools: PDFs and PDFs_lowQ grids are different'
          stop
       end if
       if (PDFs%Q_vals(i) > LHAPDF_cutoff) exit
       PDFs%tab(:,:,i) = PDFs_lowQ%tab(:,:,i)
    end do

    !!write(*,*) ' '
    !!write(*,*) ' '
    !!write(*,*) '# hybrid PDF:'
    !!do itest = 1, 100
    !!   Qtest = itest/20._dp
    !!   call EvalPdfTable_xQ(PDFs,xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!stop
    
    return
  end subroutine init_pdfs_NNLOPS

  
  !=======================================================================================
  !! initialise Hoppet's PDF table from LHAPDF
  subroutine init_pdfs_from_LHAPDF(pdf_name, pdf_set, cutoff2, asmz)
    character(len=*), intent(in)   :: pdf_name
    integer,          intent(in)   :: pdf_set
    real(dp), optional, intent(in) :: asmz, cutoff2
    !--------------------------------------------------
    real(dp) :: alphasMZ, alphasPDF, PDF_cutoff2, Qmax
    !--------------------------------------------------
    !>> for tests:
    real(dp) :: xf(-6:6), xtest, Qtest 
    integer  :: itest
    !--------------------------------------------------
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    !------------------------------------------------------------
    ! set up LHAPDF
    call InitPDFsetByName(trim(pdf_name))
    call InitPDF(pdf_set)

    ! read PDF cutoff and initialise splitting functions
    if (present(cutoff2)) then
       PDF_cutoff = sqrt(cutoff2)
    else
       call GetQ2min(pdf_set,PDF_cutoff2)
       PDF_cutoff = sqrt(PDF_cutoff2) * 2d0
    end if

    write(*,*) '==============================='
    write(*,*) 'init_pdfs_from_LHAPDF called'
    write(*,*) 'PDF_cutoff [GeV] = ',PDF_cutoff
    write(*,*) '==============================='

  
    Qmax = max(2e4_dp, cs%rts)
    call init_grid_and_dglap(Qmax)


    ! sort out the coupling: get it from LHAPDF, which will have been initialised in init_pdfs
    if (present(asmz)) then
       alphasMZ = asmz
    else
       alphasMZ = alphasPDF(MZ)
    endif

    !call InitRunningCoupling(alfas=alphasMZ, Q=MZ, nloop=3, fixnf=nf_int)
    ! Initialise with variable number of active flavours
    call InitRunningCoupling(alfas=alphasMZ, Q=MZ, nloop=3)

    call FillPdfTable_LHAPDF(PDFs, evolvePDF)

    !>> compare to LHAPDF evolution
    !!xtest = 0.5_dp
    !!write(*,*) '# test at x = ', xtest
    !!write(*,*) '# LHAPDF:'
    !!do itest = 1, 100
    !!   Qtest = itest/50._dp
    !!   call evolvePDF(xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!write(*,*) ' '
    !!write(*,*) ' '
    !!write(*,*) '# hoppet:'
    !!do itest = 1, 100
    !!   Qtest = itest/50._dp
    !!   call EvalPdfTable_xQ(PDFs,xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!stop
  end subroutine init_pdfs_from_LHAPDF


  !=======================================================================================
  !! initialise Hoppet's PDF table by evolution from an LHAPDF input
  !! at some reference scale.
  subroutine init_pdfs_5FS_evolution(pdf_name, pdf_set, cutoff2)
    character(len=*), intent(in)   :: pdf_name
    integer,          intent(in)   :: pdf_set
    real(dp), optional, intent(in) :: cutoff2
    !--------------------------------------------------
    real(dp), pointer :: pdf_Q(:,:)
    real(dp) :: Qmax, ref_Q, alphas_Q, alphasPDF, PDF_cutoff2
    !--------------------------------------------------
    !>> for tests:
    real(dp) :: xf(-6:6), xtest, Qtest 
    integer  :: itest
    !--------------------------------------------------

    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    ! set up LHAPDF
    call InitPDFsetByName(trim(pdf_name))
    call InitPDF(pdf_set)
    
    ! set PDF cutoff for hoppet grids
    if (present(cutoff2)) then
       PDF_cutoff = sqrt(cutoff2)
    else
       call GetQ2min(pdf_set,PDF_cutoff2)
       PDF_cutoff = sqrt(PDF_cutoff2) * 1.5d0
    end if
    
    ! set up the grids
    Qmax = max(2e4_dp, cs%rts)
    call init_grid_and_dglap(Qmax)
    call AllocPDF(grid, pdf_Q)

    ! extract the PDF at our reference (low) Q value
    ref_Q  = MZ            !>> do backward evolution with 5 flavours from MZ to the low scale
    alphas_Q = alphasPDF(ref_Q)
    call InitPDF_LHAPDF(grid, pdf_Q, evolvePDF, ref_Q)

    ! set up a running coupling, set up the PDF tabulation and perform the evolution
    !>> fixed flavour number scheme
    call InitRunningCoupling(alfas=alphas_Q, Q=ref_Q, nloop=3, fixnf=nf_int)

    call EvolvePdfTable(PDFs, ref_Q, pdf_Q, dglap_h, ash_global)

    !>> compare to LHAPDF evolution
    !!xtest = 0.5_dp
    !!write(*,*) '# test at x = ', xtest
    !!write(*,*) '# LHAPDF:'
    !!do itest = 1, 100
    !!   Qtest = itest/50._dp
    !!   call evolvePDF(xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!write(*,*) ' '
    !!write(*,*) ' '
    !!write(*,*) '# hoppet:'
    !!do itest = 1, 100
    !!   Qtest = itest/50._dp
    !!   call EvalPdfTable_xQ(PDFs,xtest,Qtest,xf)
    !!   write(*,*) xtest, Qtest, xf
    !!end do
    !!stop

    !------------------------------------------------------------
  end subroutine init_pdfs_5FS_evolution


  !=======================================================================================
  ! Use dummy pdfs for tests
  subroutine init_dummy_pdfs(ref_Q, alphas_Q)
    real(dp),         intent(in) :: ref_Q, alphas_Q
    !--------------------------------------------------
    real(dp), pointer :: pdf_Q(:,:)
    real(dp) :: Qmax
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    Qmax = max(2e4_dp, cs%rts)
    PDF_cutoff = 0.5_dp
    call init_grid_and_dglap(Qmax)
    call AllocPDF(grid, pdf_Q)

    pdf_Q = unpolarized_dummy_pdf(xValues(grid))

    ! set up a running coupling, set up the PDF tabulation and perform the evolution
    call InitRunningCoupling(alfas=alphas_Q, Q=ref_Q, nloop=3, fixnf=nf_int)
    ! Initialise with variable number of active flavours
    !call InitRunningCoupling(alfas=alphas_Q, Q=ref_Q, nloop=3)
    call EvolvePdfTable(PDFs, ref_Q, pdf_Q, dglap_h, ash_global)
    !------------------------------------------------------------
  end subroutine init_dummy_pdfs

  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls

    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    ! (not actually needed after setting pdf=0
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)

    ! labels iflv_g, etc., come from the hoppet_v1 module, inherited
    ! from the main program
    pdf(:, iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar

  end function unpolarized_dummy_pdf



  ! initialise coefficient functions used for NNLL luminosity
  subroutine init_C_matxs()
    implicit none
    real(dp) :: Qmax, muF_ref
    integer  :: iQ
    real(dp) :: out(-6:6)
    !----------------------------------------------
    
    call AllocSplitMat(grid, P0_matrix, nf_lcl)
    call AllocSplitMat(grid, P1_matrix, nf_lcl)
    call AllocSplitMat(grid, P2_matrix, nf_lcl)
    call InitSplitMat(P0_matrix, dglap_h%P_LO)
    call InitSplitMat(P1_matrix, dglap_h%P_NLO)
    call InitSplitMat(P2_matrix, dglap_h%P_NNLO)

    call AllocSplitMat(grid, P0_x_P0_matrix, nf_lcl)
    call AllocSplitMat(grid, C1_x_P0_matrix, nf_lcl)
    call SetToConvolution(P0_x_P0_matrix, dglap_h%P_LO, dglap_h%P_LO)
    call SetToConvolution(C1_x_P0_matrix, C1_matrix, dglap_h%P_LO) 

    call AllocSplitMat(grid, P1_x_P0_matrix, nf_lcl)
    call SetToConvolution(P1_x_P0_matrix, dglap_h%P_NLO, dglap_h%P_LO) 
    call AllocSplitMat(grid, P0_x_P1_matrix, nf_lcl)
    call SetToConvolution(P0_x_P1_matrix, dglap_h%P_LO, dglap_h%P_NLO) 
    call AllocSplitMat(grid, C1_x_P1_matrix, nf_lcl)
    call SetToConvolution(C1_x_P1_matrix,  C1_matrix, dglap_h%P_NLO) 
    
    ! O(as) coefficient functions with spin correlations don't receive any resummation
    ! scale dependence, therefore they remain unchanged

    ! In the following create a pseudo-pdf which contains the convolution between the
    ! coefficient functions and the pdfs
    ! build the PDF grid: WARNING: the following parameters are to be set as in init_grid_and_dglap
    Qmax = max(2e4_dp, cs%rts)
    muF_ref=PDF_cutoff

    call AllocPdfTable(grid, C1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, C2_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, C1_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, C2_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, G1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, P1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, P2_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, P0_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, P0_P1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, P1_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    
    call AllocPdfTable(grid, P0_P0_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, C1_P0_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, C1_P1_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, C1_P0_P0_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, C1_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, C1_P1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, C2_P0_matrix_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, C2_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)

    call AllocPdfTable(grid, G1_P0_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    call AllocPdfTable(grid, G1_P1_PDFs, muF_ref, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    
    ! build useful convolution between pdfs and coefficient functions
    do iQ=0, PDFs%nQ
       ! Start from O(as)
       C1_matrix_PDFs%tab(:,:,iQ) = C1_matrix    .conv. PDFs%tab(:,:,iQ)
       P0_PDFs%tab(:,:,iQ)        = dglap_h%P_LO .conv. PDFs%tab(:,:,iQ)
       
       C1_PDFs%tab(:,:,iQ)        = C1_matrix_PDFs%tab(:,:,iQ) + P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! C1_plus_P .conv. PDFs%tab(:,:,iQ)

       ! Now handle O(as^2)
       C2_matrix_PDFs%tab(:,:,iQ)        = C2_matrix      .conv. PDFs%tab(:,:,iQ)
       P1_PDFs%tab(:,:,iQ)               = dglap_h%P_NLO  .conv. PDFs%tab(:,:,iQ)
       P2_PDFs%tab(:,:,iQ)               = dglap_h%P_NNLO .conv. PDFs%tab(:,:,iQ)
       P0_P0_PDFs%tab(:,:,iQ)            = P0_x_P0_matrix .conv. PDFs%tab(:,:,iQ)
       C1_P0_matrix_PDFs%tab(:,:,iQ)     = C1_x_P0_matrix .conv. PDFs%tab(:,:,iQ)
       
       C2_PDFs%tab(:,:,iQ)        = C2_matrix_PDFs%tab(:,:,iQ)  &
            & + P0_PDFs%tab(:,:,iQ) * (pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2 &
            & - two*pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)*cs%ln_Q2_muR2) &
            & + P1_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + P0_P0_PDFs%tab(:,:,iQ) * (half*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2) &
            & + C1_P0_matrix_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + C1_matrix_PDFs%tab(:,:,iQ) * (-two*pi*beta0*cs%ln_Q2_muR2)  & ! C2_plus_P .conv. PDFs%tab(:,:,iQ)
            ! add the new term coming from the double emission integral (MiNNLO only)
            & + P0_PDFs%tab(:,:,iQ) * (-two*zeta3*A(1))
            
       
       G1_PDFs%tab(:,:,iQ)        = G1_matrix .conv. PDFs%tab(:,:,iQ)
       
       
       P0_P1_PDFs%tab(:,:,iQ)            = P0_x_P1_matrix .conv. PDFs%tab(:,:,iQ)
       P1_P0_PDFs%tab(:,:,iQ)            = P1_x_P0_matrix .conv. PDFs%tab(:,:,iQ)
       C1_P1_matrix_PDFs%tab(:,:,iQ)     = C1_x_P1_matrix .conv. PDFs%tab(:,:,iQ)
       
       C1_P0_PDFs%tab(:,:,iQ)     = C1_P0_matrix_PDFs%tab(:,:,iQ) + P0_P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! dglap_h%P_LO .conv. C1_PDFs%tab(:,:,iQ)
       C1_P1_PDFs%tab(:,:,iQ)     = C1_P1_matrix_PDFs%tab(:,:,iQ) + P0_P1_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! dglap_h%P_NLO .conv. C1_PDFs%tab(:,:,iQ)
       
       C2_P0_matrix_PDFs%tab(:,:,iQ)     = C2_matrix .conv. P0_PDFs%tab(:,:,iQ)
       P0_P0_P0_PDFs%tab(:,:,iQ)         = P0_matrix .conv. P0_P0_PDFs%tab(:,:,iQ)
       C1_P0_P0_matrix_PDFs%tab(:,:,iQ)  = C1_matrix .conv. P0_P0_PDFs%tab(:,:,iQ)
       !       The line below redefines C1_P1_matrix_PDFs%tab with a different convolution order
       !       C1_P1_matrix_PDFs%tab(:,:,iQ)     = C1_matrix .conv. P1_PDFs%tab(:,:,iQ)

         
       C2_P0_PDFs%tab(:,:,iQ) = C2_P0_matrix_PDFs%tab(:,:,iQ) + P0_P0_PDFs%tab(:,:,iQ) * (pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2 &
            & - two*pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)*cs%ln_Q2_muR2) &
            & + P1_P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + P0_P0_P0_PDFs%tab(:,:,iQ) * (half*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2) &
            & + C1_P0_P0_matrix_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + C1_P0_matrix_PDFs%tab(:,:,iQ) * (-two*pi*beta0*cs%ln_Q2_muR2) & ! dglap_h%P_LO .conv. C2_PDFs%tab(:,:,iQ)
            ! add the new term coming from the double emission integral (MiNNLO only)
            & + P0_P0_PDFs%tab(:,:,iQ) * (-two*zeta3*A(1))
       
       G1_P0_PDFs%tab(:,:,iQ) = G1_matrix .conv. P0_PDFs%tab(:,:,iQ)
       G1_P1_PDFs%tab(:,:,iQ) = G1_matrix .conv. P1_PDFs%tab(:,:,iQ)
    end do

    ! Useful to debug
    ! call EvalPdfTable_xQ(P0_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'P0_pdfs, init ',out
    ! call EvalPdfTable_xQ(C2_matrix_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'C2_matrix, init ',out
    ! ! >> PM: check with RadISH (add same lines in init_C_matxs of that code)
    ! call EvalPdfTable_xQ(C2_PDFs,0.5_dp,10._dp,out)  
    ! write(*,*) 'C2_pdfs,z3,a1 init ',out,zeta3,A(1)

    return
  end subroutine init_C_matxs

  ! reset coefficient functions used for NNLL luminosity
  subroutine reset_dynamical_C_matxs()
    implicit none
    integer  :: iQ
    real(dp) :: out(-6:6)
    type(process_and_parameters), save :: cs_sav
    logical, save :: ini=.true.
    !----------------------------------------------
    if(ini) then
       ini = .false.
    else
       if(compare_process_and_parameters(cs,cs_sav)) then
          return
       endif
    endif

    cs_sav = cs
    
    ! update the various grids
    do iQ=0, PDFs%nQ

       C1_PDFs%tab(:,:,iQ)        = C1_matrix_PDFs%tab(:,:,iQ) + P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! C1_plus_P .conv. PDFs%tab(:,:,iQ)

       C2_PDFs%tab(:,:,iQ)        = C2_matrix_PDFs%tab(:,:,iQ)   &
            & + P0_PDFs%tab(:,:,iQ) * (pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2 &
            & - two*pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)*cs%ln_Q2_muR2) &
            & + P1_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + P0_P0_PDFs%tab(:,:,iQ) * (half*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2) &
            & + C1_P0_matrix_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + C1_matrix_PDFs%tab(:,:,iQ) * (-two*pi*beta0*cs%ln_Q2_muR2) & ! C2_plus_P .conv. PDFs%tab(:,:,iQ)
            ! add the new term coming from the double emission integral (MiNNLO only)
            & + P0_PDFs%tab(:,:,iQ) * (-two*zeta3*A(1))


       C1_P0_PDFs%tab(:,:,iQ)     = C1_P0_matrix_PDFs%tab(:,:,iQ) + P0_P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! dglap_h%P_LO .conv. C1_PDFs%tab(:,:,iQ)
       C1_P1_PDFs%tab(:,:,iQ)     = C1_P1_matrix_PDFs%tab(:,:,iQ) + P0_P1_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) ! dglap_h%P_NLO .conv. C1_PDFs%tab(:,:,iQ)

       C2_P0_PDFs%tab(:,:,iQ) = C2_P0_matrix_PDFs%tab(:,:,iQ) + P0_P0_PDFs%tab(:,:,iQ) * (pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2 &
            & - two*pi*beta0*(-cs%ln_muF2_M2+cs%ln_Q2_M2)*cs%ln_Q2_muR2) &
            & + P1_P0_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + P0_P0_P0_PDFs%tab(:,:,iQ) * (half*(-cs%ln_muF2_M2+cs%ln_Q2_M2)**2) &
            & + C1_P0_P0_matrix_PDFs%tab(:,:,iQ) * (-cs%ln_muF2_M2+cs%ln_Q2_M2) &
            & + C1_P0_matrix_PDFs%tab(:,:,iQ) * (-two*pi*beta0*cs%ln_Q2_muR2) & ! dglap_h%P_LO .conv. C2_PDFs%tab(:,:,iQ)
            ! add the new term coming from the double emission integral (MiNNLO only)
            & + P0_P0_PDFs%tab(:,:,iQ) * (-two*zeta3*A(1))

    end do

    ! Useful to debug
    ! call EvalPdfTable_xQ(P0_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'P0_pdfs, reset ',out
    ! call EvalPdfTable_xQ(C2_matrix_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'C2_matrix, reset ',out
    ! ! >> PM: check with RadISH (add same lines in init_C_matxs of that code)
    ! call EvalPdfTable_xQ(C2_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'C2_pdfs,z3,a1 reset ',out,zeta3,A(1)

    
    return
  end subroutine reset_dynamical_C_matxs

  !======================================================================
  ! returns x*f(x) at mu2, x for POWHEG
  subroutine xfxq2_hoppet(x, mu2, xfx)
    real(dp),         intent(in)  :: x, mu2
    real(dp),         intent(out) :: xfx(-6:6)
    !---------------------------------------------------------
    real(dp) :: x_pdf(-6:7), muF
    
    x_pdf = zero
    xfx = zero
    muF = sqrt(mu2)
    
    ! evaluate pdfs at a given scale and x    
    call EvalPdfTable_xQ(PDFs, x, muF, x_pdf)

    xfx = x_pdf(-6:6)
    return
  end subroutine xfxq2_hoppet


  function alphas_hoppet(mu) result(res)
    real(dp), intent(in) :: mu
    real(dp) :: res

    if (mu > PDF_cutoff) then
       res = RunningCoupling(mu)
    else
       res = RunningCoupling(PDF_cutoff)
    end if


  end function alphas_hoppet


  
  !======================================================================
  ! return the two incoming PDFs evaluated at scale muF, and x, taking into
  ! account whether the collider is pp or ppbar
  subroutine get_pdfs_Born(muF, x1, x2, collider, pdf1, pdf2)
    real(dp),         intent(in)  :: muF, x1, x2
    character(len=*), intent(in)  :: collider
    real(dp),         intent(out) :: pdf1(-6:), pdf2(-6:)
    !---------------------------------------------------------

    pdf1=zero
    pdf2=zero
    ! evaluate pdfs at a given scale and x    
    call EvalPdfTable_xQ(PDFs, x1, muF, pdf1)
    call EvalPdfTable_xQ(PDFs, x2, muF, pdf2)
    
    select case(trim(collider))
    case("pp")
       pdf2 = pdf2
    case("ppbar")
       pdf2(-6:6) = pdf2(6:-6:-1)
       pdf2(7)    = pdf2(7)          ! index 7 contains info on representation in flavour space
    case default
       call wae_error("unrecognized collider: "//collider)
    end select
    ! Hoppet returns x*f(x)
    pdf1 = pdf1/x1
    pdf2 = pdf2/x2

  end subroutine get_pdfs_Born


  !======================================================================
  ! this works out the partonic luminosity for the LO of process proc;
  ! gives a luminosity based on the two input PDFs, which are evaluated at
  ! a given muF and x depending on the Born kinematics
  ! Remember this reads f(x) as an input instead of the conventional x*f(x) of hoppet
  function lumi_Born(pdf1, pdf2, msqB) result(res)
    real(dp),                     intent(in) :: pdf1(-6:), pdf2(-6:)
    real(dp), intent(in) :: msqB(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl) ! Matrix containing the Born matrix elements squared
    real(dp)                     :: res
    !----------------------------
    integer :: i, j

    ! build parton luminosity
    res=zero
    do i= -nf_lcl, nf_lcl
       do j= -nf_lcl, nf_lcl
          if (msqB(i,j) .ne. 0d0) then
             res = res + msqB(i,j)*pdf1(i)*pdf2(j)
          endif
       end do
    end do

  end function lumi_Born


  ! The following function returns the third order contribution D3 to the NNLOPS counterterm
  subroutine Dterms(D1, D2, D3, exp_minus_L, x1, x2, msqB, msqV1, msqV2, alphas_in, alphas_sudakov)
    real(dp),                  intent(in) :: exp_minus_L, x1, x2
    real(dp), optional, intent(in) :: alphas_in, alphas_sudakov
    real(dp), intent(in) :: msqB(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl), msqV1(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl), msqV2(-nf_lcl:nf_lcl,-nf_lcl:nf_lcl)
    real(dp),                  intent(out) :: D1, D2, D3
    !----------------------------------------------
    real(dp) :: pdf1(-6:7), pdf2(-6:7)
    real(dp) :: P0_conv_pdf1(-6:7), P0_conv_pdf2(-6:7)
    real(dp) :: P1_conv_pdf1(-6:7), P1_conv_pdf2(-6:7)
    real(dp) :: P2_conv_pdf1(-6:7), P2_conv_pdf2(-6:7)
    real(dp) :: P_conv_pdf1(-6:7), P_conv_pdf2(-6:7)
    real(dp) :: C1_conv_pdf1(-6:7), C1_conv_pdf2(-6:7)
    real(dp) :: C2_conv_pdf1(-6:7), C2_conv_pdf2(-6:7)
    real(dp) :: G1_conv_pdf1(-6:7), G1_conv_pdf2(-6:7)
    real(dp) :: P0_C1_conv_pdf1(-6:7), P0_C1_conv_pdf2(-6:7)
    real(dp) :: P0_C2_conv_pdf1(-6:7), P0_C2_conv_pdf2(-6:7)
    real(dp) :: P0_G1_conv_pdf1(-6:7), P0_G1_conv_pdf2(-6:7)
    real(dp) :: P1_C1_conv_pdf1(-6:7), P1_C1_conv_pdf2(-6:7)
    real(dp) :: P1_G1_conv_pdf1(-6:7), P1_G1_conv_pdf2(-6:7)
    real(dp) :: P_C1_conv_pdf1(-6:7), P_C1_conv_pdf2(-6:7)
    real(dp) :: P_C2_conv_pdf1(-6:7), P_C2_conv_pdf2(-6:7)
    real(dp) :: P_G1_conv_pdf1(-6:7), P_G1_conv_pdf2(-6:7)
    integer  :: i
    real(dp) :: as2pi, as2pi_sudakov, L, lambda, as2pi_M, muF, muR, H1
    real(dp) :: dS1, dS2, dS3, L0, L1, L2, dL1, dL2, dL3 ! building blocks of D3
    !as in appendix C of the paper >> the following flag must be set
    !to false is using init_pdfs_NNLOPS() for consistency with POWHEG
    logical, parameter :: use_analytic_alpha = .false. 
    real(dp) :: out(-6:6)


    D1 = zero; D2 = zero; D3 = zero

    muF = cs%muF * exp_minus_L ! when pt << M this becomes cs%muF * exp(-log(Q/pT)) = pT * cs%muF/Q = pT * KF
    muR = cs%muR * exp_minus_L
    ! define logarithm for expansion of the Sudakov
    L = - log(exp_minus_L)

    ! get the pdf grids and the various convolutions with the coefficient functions
    P0_conv_pdf1=zero
    P0_conv_pdf2=zero
    call EvalPdfTable_xQ(P0_PDFs, x1, muF, P0_conv_pdf1)
    call EvalPdfTable_xQ(P0_PDFs, x2, muF, P0_conv_pdf2)
    ! Hoppet returns x*f(x)
    P0_conv_pdf1 = P0_conv_pdf1/x1
    P0_conv_pdf2 = P0_conv_pdf2/x2
    ! ---
    P1_conv_pdf1=zero
    P1_conv_pdf2=zero
    call EvalPdfTable_xQ(P1_PDFs, x1, muF, P1_conv_pdf1)
    call EvalPdfTable_xQ(P1_PDFs, x2, muF, P1_conv_pdf2)
    P1_conv_pdf1 = P1_conv_pdf1/x1
    P1_conv_pdf2 = P1_conv_pdf2/x2
    ! ---
    P2_conv_pdf1=zero
    P2_conv_pdf2=zero
    call EvalPdfTable_xQ(P2_PDFs, x1, muF, P2_conv_pdf1)
    call EvalPdfTable_xQ(P2_PDFs, x2, muF, P2_conv_pdf2)
    P2_conv_pdf1 = P2_conv_pdf1/x1
    P2_conv_pdf2 = P2_conv_pdf2/x2

    ! ---
    C1_conv_pdf1=zero
    C1_conv_pdf2=zero
    call EvalPdfTable_xQ(C1_PDFs, x1, muF, C1_conv_pdf1)
    call EvalPdfTable_xQ(C1_PDFs, x2, muF, C1_conv_pdf2)
    C1_conv_pdf1 = C1_conv_pdf1/x1
    C1_conv_pdf2 = C1_conv_pdf2/x2
    ! ---
    C2_conv_pdf1=zero
    C2_conv_pdf2=zero
    call EvalPdfTable_xQ(C2_PDFs, x1, muF, C2_conv_pdf1)
    call EvalPdfTable_xQ(C2_PDFs, x2, muF, C2_conv_pdf2)
    C2_conv_pdf1 = C2_conv_pdf1/x1
    C2_conv_pdf2 = C2_conv_pdf2/x2
    ! ---
    G1_conv_pdf1=zero
    G1_conv_pdf2=zero
    call EvalPdfTable_xQ(G1_PDFs, x1, muF, G1_conv_pdf1)
    call EvalPdfTable_xQ(G1_PDFs, x2, muF, G1_conv_pdf2)
    G1_conv_pdf1 = G1_conv_pdf1/x1
    G1_conv_pdf2 = G1_conv_pdf2/x2

    ! ---
    P0_C1_conv_pdf1=zero
    P0_C1_conv_pdf2=zero
    call EvalPdfTable_xQ(C1_P0_PDFs, x1, muF, P0_C1_conv_pdf1)
    call EvalPdfTable_xQ(C1_P0_PDFs, x2, muF, P0_C1_conv_pdf2)
    P0_C1_conv_pdf1 = P0_C1_conv_pdf1/x1
    P0_C1_conv_pdf2 = P0_C1_conv_pdf2/x2
    ! ---
    P0_C2_conv_pdf1=zero
    P0_C2_conv_pdf2=zero
    call EvalPdfTable_xQ(C2_P0_PDFs, x1, muF, P0_C2_conv_pdf1)
    call EvalPdfTable_xQ(C2_P0_PDFs, x2, muF, P0_C2_conv_pdf2)
    P0_C2_conv_pdf1 = P0_C2_conv_pdf1/x1
    P0_C2_conv_pdf2 = P0_C2_conv_pdf2/x2
    ! ---
    P0_G1_conv_pdf1=zero
    P0_G1_conv_pdf2=zero
    call EvalPdfTable_xQ(G1_P0_PDFs, x1, muF, P0_G1_conv_pdf1)
    call EvalPdfTable_xQ(G1_P0_PDFs, x2, muF, P0_G1_conv_pdf2)
    P0_G1_conv_pdf1 = P0_G1_conv_pdf1/x1
    P0_G1_conv_pdf2 = P0_G1_conv_pdf2/x2


    ! ---
    P1_C1_conv_pdf1=zero
    P1_C1_conv_pdf2=zero
    call EvalPdfTable_xQ(C1_P1_PDFs, x1, muF, P1_C1_conv_pdf1)
    call EvalPdfTable_xQ(C1_P1_PDFs, x2, muF, P1_C1_conv_pdf2)
    P1_C1_conv_pdf1 = P1_C1_conv_pdf1/x1
    P1_C1_conv_pdf2 = P1_C1_conv_pdf2/x2
    ! ---
    P1_G1_conv_pdf1=zero
    P1_G1_conv_pdf2=zero
    call EvalPdfTable_xQ(G1_P1_PDFs, x1, muF, P1_G1_conv_pdf1)
    call EvalPdfTable_xQ(G1_P1_PDFs, x2, muF, P1_G1_conv_pdf2)
    P1_G1_conv_pdf1 = P1_G1_conv_pdf1/x1
    P1_G1_conv_pdf2 = P1_G1_conv_pdf2/x2

    ! Define running coupling
    if(present(alphas_in)) then
       as2pi = alphas_in/twopi
    else
       if (.not.use_analytic_alpha) then
          if (muR > PDF_cutoff) then
             as2pi = RunningCoupling(muR)/twopi
          else
             as2pi = RunningCoupling(PDF_cutoff)/twopi
          end if
       else
          lambda = RunningCoupling(cs%M)*beta0*L
          if (lambda < half) then
             as2pi = RunningCoupling(cs%M)/twopi
             as2pi = as2pi / (1-two*lambda) * (one &
                  & - (twopi*as2pi) / (1-two*lambda) * beta1/beta0 * log(one-two*lambda))
          else
             write(*,*) 'WARNING: freezing alphas at the Landau singularity'
             as2pi = cs%alphas_muR/twopi
          end if
       end if
    end if

    ! Finally, define running coupling used in the Sudakov
    as2pi_sudakov = as2pi

    
    ! Now proceed with the derivative of the luminosity
    call get_pdfs_Born(muF, x1, x2, cs%collider, pdf1, pdf2)
            

    ! create macros for various convolutions as coefficients of as2pi^i
    P0_conv_pdf1 = two*(P0_conv_pdf1)
    P0_conv_pdf2 = two*(P0_conv_pdf2)
    P1_conv_pdf1 = two*(P1_conv_pdf1)
    P1_conv_pdf2 = two*(P1_conv_pdf2)
    P2_conv_pdf1 = two*(P2_conv_pdf1)
    P2_conv_pdf2 = two*(P2_conv_pdf2)

    P0_C1_conv_pdf1 = two*(P0_C1_conv_pdf1)
    P0_C1_conv_pdf2 = two*(P0_C1_conv_pdf2)
    P1_C1_conv_pdf1 = two*(P1_C1_conv_pdf1)
    P1_C1_conv_pdf2 = two*(P1_C1_conv_pdf2)

    P0_C2_conv_pdf1 = two*(P0_C2_conv_pdf1)
    P0_C2_conv_pdf2 = two*(P0_C2_conv_pdf2)
    P0_G1_conv_pdf1 = two*(P0_G1_conv_pdf1)
    P0_G1_conv_pdf2 = two*(P0_G1_conv_pdf2)

    ! define building blocks
    dS1 = two * (two*A(1)*L + (B(1) - A(1)*cs%ln_Q2_M2)) * as2pi_sudakov
    dS2 = two * (two*A(2)*L + (B(2) - A(2)*cs%ln_Q2_M2)) * as2pi_sudakov**2
    dS3 = two * (two*A(3)*L + (B(3) - A(3)*cs%ln_Q2_M2)) * as2pi_sudakov**3

    L0  = lumi_Born(pdf1, pdf2, msqB)
    L1  = (lumi_Born(pdf1, pdf2, msqV1) + lumi_Born(C1_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, C1_conv_pdf2, msqB)) * as2pi
    L2  = (lumi_Born(pdf1, pdf2, msqV2) + lumi_Born(C2_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, C2_conv_pdf2, msqB) &
         & + lumi_Born(C1_conv_pdf1, pdf2, msqV1) + lumi_Born(pdf1, C1_conv_pdf2, msqV1) &
         & + lumi_Born(C1_conv_pdf1, C1_conv_pdf2, msqB) + lumi_Born(G1_conv_pdf1, G1_conv_pdf2, msqB)) * as2pi**2   

    dL1 = (lumi_Born(pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, pdf2, msqB)) * as2pi
    dL2 = (lumi_Born(pdf1, P1_conv_pdf2, msqB) + lumi_Born(P1_conv_pdf1, pdf2, msqB) &
         & + (lumi_Born(pdf1, P0_conv_pdf2, msqV1) + lumi_Born(P0_conv_pdf1, pdf2, msqV1)) &
         & + (lumi_Born(C1_conv_pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, C1_conv_pdf2, msqB)) &
         & + (lumi_Born(P0_C1_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, P0_C1_conv_pdf2, msqB)) &
         & - (4._dp*pi*beta0) * (lumi_Born(C1_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, C1_conv_pdf2, msqB)) &
         & - (4._dp*pi*beta0) * lumi_Born(pdf1, pdf2, msqV1)) * as2pi**2
    dL3 = (lumi_Born(P0_conv_pdf1, pdf2, msqV2)  + lumi_Born(pdf1, P0_conv_pdf2, msqV2) &
         & + lumi_Born(P1_conv_pdf1, pdf2, msqV1)  + lumi_Born(pdf1, P1_conv_pdf2, msqV1) &
         & + lumi_Born(C1_conv_pdf1, P0_conv_pdf2, msqV1) + lumi_Born(P0_conv_pdf1, C1_conv_pdf2, msqV1) &
         & + lumi_Born(pdf1, P0_C1_conv_pdf2, msqV1) + lumi_Born(P0_C1_conv_pdf1, pdf2, msqV1) &
         & + lumi_Born(P2_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, P2_conv_pdf2, msqB) &
         & + lumi_Born(C2_conv_pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, C2_conv_pdf2, msqB) &
         & + lumi_Born(P0_C2_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, P0_C2_conv_pdf2, msqB) &
         & + lumi_Born(C1_conv_pdf1, P1_conv_pdf2, msqB) + lumi_Born(P1_conv_pdf1, C1_conv_pdf2, msqB) &
         & + lumi_Born(P1_C1_conv_pdf1, pdf2, msqB) + lumi_Born(pdf1, P1_C1_conv_pdf2, msqB) &
         & + lumi_Born(P0_C1_conv_pdf1, C1_conv_pdf2, msqB) + lumi_Born(C1_conv_pdf1, P0_C1_conv_pdf2, msqB) &
         & + lumi_Born(P0_G1_conv_pdf1, G1_conv_pdf2, msqB) + lumi_Born(G1_conv_pdf1, P0_G1_conv_pdf2, msqB)) * as2pi**3 &
         & - four*twopi*beta0 * L2 * as2pi &
         & - two*twopi**2*beta1 * L1 * as2pi**2

    ! ---------------------------------- as[pt]^1 (beginning) ---------------------------------------------
    D1 = dS1*L0 + dL1
    ! ------------------------------------- as[pt]^1 (end) ------------------------------------------------


    ! ---------------------------------- as[pt]^2 (beginning) ---------------------------------------------
    D2 = dS2*L0 + dS1*L1 + dL2
    ! now add the explicit scale dependence (see NNLOPS.nb for a derivation of the formulae)
    D2 = D2 - as2pi * twopi*beta0 * dL1 * (cs%ln_muF2_M2 - cs%ln_muR2_M2)   
    ! ------------------------------------- as[pt]^2 (end) ------------------------------------------------


    ! ---------------------------- as[pt]^3 (beginning) - App. C of the paper -----------------------------
    D3 = (dS1 * L2 + dS2 * L1 + dS3 * L0 + dL3)


    ! now add the explicit scale dependence (see NNLOPS.nb for a derivation of the formulae)
    D3 = D3 + (- twopi*(twopi*beta1 * (lumi_Born(pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, pdf2, msqB)) &
         &  + beta0*(lumi_Born(pdf1, P0_conv_pdf2, msqV1) + lumi_Born(P0_conv_pdf1, pdf2, msqV1) &
         &  + two*(lumi_Born(pdf1, P1_conv_pdf2, msqB) + lumi_Born(P1_conv_pdf1, pdf2, msqB)) &
         &  + (lumi_Born(C1_conv_pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, C1_conv_pdf2, msqB)) &
         &  + (lumi_Born(pdf1, P0_C1_conv_pdf2, msqB) + lumi_Born(P0_C1_conv_pdf1, pdf2, msqB)))) * (cs%ln_muF2_M2 - cs%ln_muR2_M2) &
         &  + twopi**2*beta0**2*(lumi_Born(pdf1, P0_conv_pdf2, msqB) + lumi_Born(P0_conv_pdf1, pdf2, msqB)) &
         &  * (cs%ln_muF2_M2 - cs%ln_muR2_M2)**2) * as2pi**3
    ! ------------------------------------- as[pt]^3 (end) ------------------------------------------------

    ! Useful to debug
    ! call EvalPdfTable_xQ(P0_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'P0_pdfs, dterms ',out
    ! call EvalPdfTable_xQ(C2_matrix_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'C2_matrix, dterms ',out
    ! ! >> PM: check with RadISH (add same lines in init_C_matxs of that code)
    ! call EvalPdfTable_xQ(C2_PDFs,0.5_dp,10._dp,out)
    ! write(*,*) 'C2_pdfs,z3,a1 dterms ',out,zeta3,A(1)

    return
  end subroutine Dterms  

end module pdfs_tools
