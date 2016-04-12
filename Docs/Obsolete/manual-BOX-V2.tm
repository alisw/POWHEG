<TeXmacs|1.0.7.14>

<style|generic>

<\body>
  <doc-data|<doc-title|The <with|font-family|tt|V2> extensions of the
  <with|font-family|tt|POWHEG BOX>>>

  <section|Radiation from resonances>

  In the <with|font-family|tt|POWHEG BOX V2> it is possible to include strong
  radiation generated from decaying resonances. In order to do so, the
  <with|font-family|tt|init_processes.f> file, where the lists of Born and
  Real processes are defined, should also define the following arrays:

  <with|font-family|tt|integer flst_bornres(1:nlegborn,flst_nborn)>

  whose entries are zero for the initial state particles and for the
  particles directly produced in the reaction, while, for particles coming
  from the decay of a resonance, it is equal to the position of the resonance
  in the particle list. Thus, for example, in <math|t<wide|t|\<bar\>>>
  production and decay:

  <with|font-family|tt|flst_born(1:nlegborn,j) \ \ =[0, \ 0, \ 6, -6,
  24,-24,-11, 12, 13,-14, \ 5, \ -5]>

  <with|font-family|tt|flst_bornres(1:nlegborn,j)=[0, \ 0, \ 0, \ 0, \ 3,
  \ 4, \ 6, \ 6, \ 7, \ 7, \ 3, \ \ 4]>

  with the <math|t<wide|t|\<bar\>>> pair coming from the production process,
  with the corresponding <with|font-family|tt|bornres> entry equal to 0. The
  <math|W<rsup|+>> and the <math|b> come from the decay of the top, which is
  the 3<math|<rsup|rd>> entry, so that the bornres entry is 3; the
  <math|e<rsup|+>> and the <math|\<nu\><rsub|e>> come from the
  <math|W<rsup|+>>, which is the 6<math|<rsup|th>> entry, and so on.

  An analogous array <with|font-family|tt|flst_realres> is supplied for the
  real emission graph. In this case, the radiated parton can either arise
  from production (<with|font-family|tt|realres> entry equal to zero) or from
  a resonance (<with|font-family|tt|realres> entry equal to the resonance
  position).

  The <with|font-family|tt|BOX>, on the basis of the provided lists, sets up
  the parameter <with|font-family|tt|flst_nreson> equal to the number of
  resonances that can radiate. This includes radiation from production, that
  is treated as a fictitious resonance, indexed by 0. Furthermore, the array
  <with|font-family|tt|flst_reslist> is set up, that contains the index of
  each resonance that can radiate. Its first entry is always 0, corresponding
  to radiation in production. If there are no radiating resonances, the BOX
  sets <with|font-family|tt|flst_nreson> to 1,
  <with|font-family|tt|flst_reslist(1)=0>.

  The user must provide a <with|font-family|tt|setreal> routine that returns
  the matrix element for radiation of a specific resonance. In order to do
  so, the <with|font-family|tt|setreal> routine must inspect the variable
  <with|font-family|tt|kn_resemitter>, which is a pointer to the resonance
  that is radiating, and supply the corresponding real matrix element. The
  <with|font-family|tt|setvirtual> routine should provide the sum of the
  virtual contributions for the production and decay processes.

  <section|Electroweak processes>

  The <with|font-family|tt|POWHEG BOX> is now capable to handle certain
  electroweak processes, where soft-photon radiation is generated with the
  shower method. The following variables are relevant for this extension:

  <with|font-family|tt|flg_with_em> must be set to true in the initialization
  of the user process

  <with|font-family|tt|pdf_nparton> must be set to 22, if incoming photons
  from parton densities are considered

  <with|font-family|tt|em_muren2> \ \ \ must be set to the (squared)
  renormalization scale for the EW corrections

  <with|font-family|tt|em_alpha> \ \ \ \ must be set to
  <math|\<alpha\><rsub|em>>.

  The <with|font-family|tt|BOX> is capable of handling photon initial and
  final state radiation from charged particles. It has been presently tested
  only with <math|W> and <math|Z> production with leptonic decays.

  <\bibliography|bib|JHEP|paper.bib>
    <\bib-list|1>
      <bibitem*|1><label|bib-Melia:2011gk>T.<nbsp>Melia, P.<nbsp>Nason,
      R.<nbsp>Rontsch, and G.<nbsp>Zanderighi,
      <with|font-shape|italic|W<rsup|+>W<rsup|+> plus dijet production in the
      POWHEGBOX>, <hlink|<with|font-family|tt|1102.4846>|http://xxx.lanl.gov/abs/1102.4846>.
      * Temporary entry *.

      <bibitem*|2><label|bib-Melia:2010bm>T.<nbsp>Melia, K.<nbsp>Melnikov,
      R.<nbsp>Rontsch, and G.<nbsp>Zanderighi,
      <with|font-shape|italic|Next-to-leading order QCD predictions for
      <math|W<rsup|+>W<rsup|+>j*j> production at the LHC>,
      <with|font-shape|italic|JHEP> <with|font-series|bold|1012> (2010) 053,
      [<hlink|<with|font-family|tt|1007.5313>|http://xxx.lanl.gov/abs/1007.5313>].

      <bibitem*|3><label|bib-Alioli:2010xd>S.<nbsp>Alioli, P.<nbsp>Nason,
      C.<nbsp>Oleari, and E.<nbsp>Re, <with|font-shape|italic|A general
      framework for implementing NLO calculations in shower Monte Carlo
      programs: the POWHEG BOX>, <with|font-shape|italic|JHEP>
      <with|font-series|bold|1006> (2010) 043,
      [<hlink|<with|font-family|tt|1002.2581>|http://xxx.lanl.gov/abs/1002.2581>].
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|par-hyphen|normal>
    <associate|sfactor|5>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|1.3|3>>
    <associate|auto-5|<tuple|1.4|3>>
    <associate|auto-6|<tuple|2|3>>
    <associate|auto-7|<tuple|2|3>>
    <associate|bib-Alioli:2010xd|<tuple|3|?>>
    <associate|bib-Frixione:2007vw|<tuple|4|?>>
    <associate|bib-Melia:2010bm|<tuple|2|?>>
    <associate|bib-Melia:2011gk|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Radiation
      from resonances> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>