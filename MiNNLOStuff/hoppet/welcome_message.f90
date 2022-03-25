! keep this outside a module so that we can recompile (e.g. for
! changing version number) without modifying any dependences.
!
subroutine HoppetWelcomeMessage
  write(*,'(a)') '-----------------------------------------------------------'
  write(*,'(a)') '               Welcome to HOPPET v. 1.2.0                  '
  write(*,'(a)') '   Higher Order Perturbative Parton Evolution Toolkit      '
  write(*,'(a)') ''
  write(*,'(a)') '        Written by Gavin P. Salam (2001-2012)'
  write(*,'(a)') '          with contributions from Juan Rojo'
  write(*,'(a)') '        Frederic Dreyer and Alexander Karlberg'
  write(*,'(a)') ''
  write(*,'(a)') ' It is made available under the GNU public license,'
  write(*,'(a)') ' with the additional request that if you use it or any'
  write(*,'(a)') ' derivative of it in scientific work then you should cite:'
  write(*,'(a)') ' G.P. Salam & J. Rojo, CPC 180(2009)120 (arXiv:0804.3755).'
  write(*,'(a)') ' '
  write(*,'(a)') ' You are also encouraged to cite the original references,'
  write(*,'(a)') ' for LO, NLO and NNLO splitting functions, the QCD'
  write(*,'(a)') ' 1, 2 and 3 loop beta functions and the coupling and '
  write(*,'(a)') ' PDF and coupling mass threshold matching functions.'
  write(*,'(a)') '-----------------------------------------------------------'
end subroutine HoppetWelcomeMessage
  
