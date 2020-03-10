<TeXmacs|1.99.7>

<style|generic>

<\body>
  <section|Corrections for upper bound violations>

  \;

  In a hit-and-miss generator, suppose that the upper bound is not working.
  The usual hit-and-miss procedure will then generate some events above the
  bound. These events will be generated as if the real distribution was
  truncated at the hight of the upper bound. In order to remedy to this, the
  weight of ub violating event could be increased by a factor equal to the
  value of the function at that point divided by the upper bound. The
  algorithm becomes:

  <\itemize-dot>
    <item>Generate a random phase space point. If the value of the function
    is below the upper bound, accept/reject.

    <item>If it is in a point where the bound is violated, accept the event
    but increase its weight by a factor equal to the ratio of the value of
    the function over the upper bound <math|f/u>.
  </itemize-dot>

  \;

  <with|gr-mode|<tuple|edit|text-at>|gr-frame|<tuple|scale|1cm|<tuple|0.100038gw|0.180031gh>>|gr-geometry|<tuple|geometry|0.533352par|0.466672par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-dash-style|10|<graphics||<spline|<point|0|1.2>|<point|1.0|2.6>|<point|2.4|3.4>|<point|3.2|3.0>|<point|4.0|2.3>|<point|4.8|3.0>|<point|5.6|3.2>|<point|7.0|2.6>>|<line|<point|0|3>|<point|7.0|3.0>>|<line|<point|0|0>|<point|7.0|0.0>>|<point|0.2|2.5>|<point|0.4|1>|<point|2.4|2.6>|<text-at|miss|<point|0.2|2.6>>|<text-at|hit|<point|0.4|1.2>>|<with|dash-style|10|<line|<point|1.4|2.99996>|<point|1.4|0.0>>>|<with|dash-style|10|<line|<point|3.2|3>|<point|3.2|0.0>>>|<with|dash-style|10|<line|<point|4.8|3>|<point|4.8|0.0>>>|<with|dash-style|10|<line|<point|6.1|3.00679>|<point|6.1|0.0>>>|<text-at|A|<point|0.7|0.5>>|<text-at|B|<point|3.9|0.5>>|<text-at|C|<point|6.5|0.5>>|<text-at|hit
  with|<point|1.6|2.3>>|<text-at|ub violation|<point|1.6|2.0>>|<text-at|D|<point|2.2|1.4>>|<text-at|G|<point|3.8|2.6>>|<text-at|H|<point|6.7|2.74082592106522>>|<text-at|F|<point|0.3|2.2>>|<text-at|E|<point|5.3|1.6>>>>

  \;

  This way , if <math|N> events are generate with a cross section <math|W>,
  histogramming all events will lead to a cross section

  <\equation*>
    <frac|<around*|(|N-N<rsub|\<gtr\>>|)>W+<big|sum><rsub|N<rsub|\<gtr\>>>W
    <frac|f<rsub|i>|u><rsub|>|N>.
  </equation*>

  This factor is reported in the counter file as 'Weight increment factor due
  to corrections for upper bound violation'.

  In order to correct for this problem, one should divide the cross section
  by that factor.

  \;

  <section|Corrections for wrong estimate of the cross section.>

  If the estimate of the accumulated cross section (in absolute value!) is
  off, we get an error that will not be eliminated no matter how many events
  one generates. POWHEG provides an independent calculation of the absolute
  value of the cross section during event generation. This feature is
  independent upon the ubexcess_correct one.

  The python script FindReweightFromCounters.py gets these results out of the
  counters.

  We have two possible corrections: we may correct for the UB violations, in
  which case we divide by the 'Weight increment ...' number given by the
  script; and we may correct for the ratio of the cross section computed on
  the fly over the stage 4 cross section, if the first one is more accurate.

  The script also spits out a total correction factor. In reality, while
  there is no reason to apply the UB correction, it makes sense to also
  correct for the total cross section only if the error of the stage 4 cross
  section estimate is smaller than the cross section used for generation.

  <section|Implementation in the POWHEG integrator.>

  The implementation of the corrections mentioned above is not totally
  straightforward, in view of the several different options under which the
  integrator may work. In all cases, each coodinate is divided up into 50
  steps, and the integration hypercube is divided into
  <math|50<rsup|n<rsub|dim>>> cells. Each cell is generated with a given
  probability, and is assigned a weight equal to the inverse of the
  generation probability. In the most sofisticated one, the integrator
  determines two sets of upper bound on the integrand, one called
  <with|font|Tlwg|ymax>, that has the form of a product of step functions for
  each coordinates, that we now call <math|Y<rsub|max><around*|(|<wide|x|\<vect\>>|)>>,
  where <math|<wide|x|\<vect\>>> is a multidimensional vector in the
  integration hypercube, and the other has the form
  <math|Y<rsub|rat><around*|(|<wide|x|\<vect\>>|)>
  B<around*|(|<wide|x|\<vect\>>|)>>, where
  <math|B<around*|(|<wide|x|\<vect\>>|)>> is the Born cross section, and
  <math|Y<rsub|rat>> is again a product of step function in each coordinate.

  \;
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Corrections
      for upper bound violations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Corrections
      for wrong estimate of the cross section.>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>