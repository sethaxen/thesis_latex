Sunday:
- read ensemble papers
  - pay attention to alternative time complexities and note cases where ensemble appears actually quite broad
    - hummerBayesianEnsembleRefinement2015: both reweighting and replica-averaged restraint (todo: it looks like they don't deconvolve ensemble width and posterior width? verify that they return a single ensemble)
    - bonomiMetainferenceBayesianInference2016: metainference, it's a Bayesian replica-averaging technique. test with NMR/RDC on ubiquitin bonomiMetadynamicMetainferenceEnhanced2016 adds metadynamics for enhanced sampling. (todo: how many replicas can/do they sample? is it weighted? how do they test convergence? check supplements. How did they average NOE data? Can we use RDC's???)
    - molnarCysScanningDisulfideCrosslinking2014 PhoQ paper, it's Bayesian two-state modeling
- write discussion
- finalize benchmark
- write final paragraph of intro
- write abstract
- address notes
Monday:
- writing session with Andrej
- finalize title
- resolve Andrej's notes
- finish appendix
Tuesday
- send Chapter 3/Appendix A to commitee
- write SHG chapter/appendix
- write long-form thesis abstract
Wednesday
- send SHG chapter/appendix and intro to commitee



https://docs.google.com/document/d/1yvZOU5Wm1guKBYY_veWHWLrpM20DVW_Zt5sRxDy6C_k/edit

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2935091/pdf/nihms225720.pdf
    has a discussion at the end about ways to apply the harmonic approach to rigid body systems in proteins

Need refs:
    Different types of "ensemble modeling" (check IMP proposal):
        multi-state modeling
        temporal ensemble? (using an MD trajectory or coupled trajectories)
        ensemble reweighting
        continuum representations?
            https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00756
                smooth continuous dielectric function, not structural modeling
        others???
        https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1006641
            Trewhella paper
    
    NOE ensemble average
        average over ensemble (spatial ensemble, not just time)
            - origin of <r^{-6}>
    
    Order parameters??
        Simultaneous determination of protein structure and dynamics (NOEs and order parameters over ensembles)

    Bayesian inference
    probabilistic programming on proteins:
    https://ieeexplore.ieee.org/document/8791469
        Probabilistic programming in Pyro for protein structure superposition
    https://pubs.acs.org/doi/10.1021/acs.jpca.0c05026
        Use PyMC3 to fit parameters of DEER distribution. No structural DOFs?
    https://www-jstor-org.ucsf.idm.oclc.org/stable/pdf/25053290.pdf
        protein structure prediction as probabilistic programming


    Misc:
        Protein structure validation and refinement using amide proton chemical shifts derived from quantum mechanics. (uses HMC?)

Flow:
    We need to fit data derived from ensembles.
    In particular, from multimodal distributions with variation about the modes.

    We want to propagate uncertainty in the data and prior information to the parameters of the ensemble. (inference)
    Probabilistic programming languages enable inference.

    There are existing approaches.
    They are either not appropriate for inference or for fitting continuum models.

    We develop a "kinematic ensemble representation" that can be integrated into existing probabilistic programming pipelines, including the integrative modeling platform.


questions re review:


do they fit a single ensemble?
do they return an ensemble of ensembles?
or do they return a sample of structures from an ensemble of ensembles?

Pros of existing approaches:
    Operate on atomic and coarse-grain representations, leveraging all the machinery of MD and enhanced sampling techniques.

You can fit a single structure, unless your system is very hetereogenous and your data reflects that at the structural resolution you are modeling.
If your distribution is a mixture of narrow distributions, you can fit multiple replicas.
If it's too inefficient to sample, you can draw a sample of structures using some technique and then assign weights to the structures, mixture weights, gives a single ensemble, only works well if the system's does not deviate very far from the initial sample (really, should be contained in).
If there are many modes or if the components are wide, then you can add many components. This scales as O(n!) to sample unless you discard indices, in which case you're not separating uncertainty and heterogeneity.
The true distribution is high-dimensional and highly multimodal. However, for many proteins, there will be large peaks containing many sub-peaks. The multi-state approach essentially approximates this distribution locally at fixed points. On the upside, it can capture sharp information like the exact location of a well, but it is not able to capture more smooth information, it needs many replicas.
This can be though of as analogous to a zero-order Taylor series approximation of the density function of the ensemble, where the function is approximated only at fixed points.
We take a different approach to approximating the same density.
Our approach approximates first the low-resolution features of the smooth ensemble; as the complexity increases, higher resolution features may be approximated.

In particular, large complex systems are typically not amenable to modeling at high resolution.
In many cases, a low resolution model will be sufficient. However, we do not always have a low resolution forward model.
Our approach permits one to maintain an atomic resolution structure but reduce the degrees of freedom using rigid bodies. The forward model can still be computed at the atomic resolution, but by taking expectations over the smooth ensemble.

data can be sensitive to low population levels, i.e. to the tails of the ensemble distribution.
A replica averaged approach could require a very large number of replicas to satisfy such data. However, in our approach, tails are explicitly modeled, hence data restraints can influence the ensemble early. However, tails are still subject to rigid body approximation.

this gives us coordinates wrt one body. 



Hi everyone,

    Attached is a new draft that primarily adds to the appendix of the SHG/TPF section. Right now the appendix lays out a rough draft of the content I would like to cover. All theory necessary to use SHG/TPF with the kinematic ensemble approach is there, it just needs to flow better, and there are a few notes for explaining the motivation. In addition to the general arguments in the NOE chapter, for SHG/TPF, there's an excellent case to be made from both the information content of the data and from the perspective of computational efficiency for using the kinematic ensemble approach.

    I'll be adding text giving the context of what SHG/TPF has been used for in structural biology, what it uniquely gives us, and what it doesn't give us. Right now this is all interleaved in the appendix, but the rationale and context will be split out to form Chapter 4, with the math remaining in the appendix.
    
    Seth

    For my records, this is git commit e2a9f6c
