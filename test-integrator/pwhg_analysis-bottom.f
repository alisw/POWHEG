c  The next subroutines, open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones
c  init   :  opens the histograms
c  topout :  closes them
c  pwhgfill  :  fills the histograms with data

      subroutine init_hist
      implicit none


      call inihists

      call bookupeqbins('tot-abs',1d0,-0.5d0,0.5d0)
      call bookupeqbins('tot',1d0,-0.5d0,0.5d0)

      end

      subroutine analysis(dsig0)
      implicit none
      real * 8  dsig0,dsig(10)

      call multi_plot_setup(dsig0,dsig,10)
      call filld('tot',0d0,dsig)
      call filld('tot-abs',0d0,abs(dsig))

      end
