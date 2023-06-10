Usage
=====

.. _paranome:

Whole paranome delineation
------------

The first and uppermost step of constructing *K*\ :sub:`S` distribution is to infer the paralogous gene families (whole paranome). The simplest way of achieving it is as below.

.. code-block:: console

   (ENV) $ wgd dmd cds.fasta

This step is actually affected by many parameters. As shown below, two key parameters are the ``eval``, which determines the e-value cut-off for similarity in diamond and thus directly impacts the whole query-subject hits table, and the ``inflation``, which determines the inflation factor in the MCL program and thus directly affects the gene family clustering result. Some other minor (maybe also major in some conditions) parameters including the ``normalizedpercent`` controling the percentage of upper hits used for gene length normalization and the ``bins`` controlling the number of bins divided in gene length normalization, both of which might influence the normalization result and thus the gene family clustering.

.. _mrbh:

Global/Local MRBHs inference
------------

RBHs are good candidates for representing orthologous relationship between genes. The global MRBHs (Multiple RBHs) can be inferred with the code below.

.. code-block:: console

   (ENV) $ wgd dmd cds1.fasta cds2.fasta cds3.fasta --globalmrbh

The parameter that can directly affect the result of global MRBHs is the ``cscore``. If it's set as default None, the process will go as strict RBHs. If it's set as a float between 0 and 1, the number of inter-specific gene pairs will probably increase, so does the number of global MRBHs. Normally, more species lead to less global MRBHs. Sometimes lowering the ``cscore`` might be needed to acquire enough global MRBHs for other analysis such as phylogenetic inference.

The local MRBHs which focus on a specific species, is a logical regime for searching the pair-wise orthologues to build orthogroups used in phylogenetic dating. The following command can be used to infer local MRBHs. The influential parameters are the same as global MRBHs.

.. code-block:: console

   (ENV) $ wgd dmd cds1.fasta cds2.fasta cds3.fasta --focus cds1.fasta

.. _orthogroup:

Orthogroup inference
------------

Orthogroups are the most common units used in evolutionary biology for comparative analysis. It can be inferred by the command below.

.. code-block:: console

   (ENV) $ wgd dmd cds1.fasta cds2.fasta cds3.fasta -oi

There are two ways of inferring orthogroups: one is the concatenation method which firstly concatenates all cds files as a super cds file and secondly infers the query-subject hits table as the single genome does, with gene length bias correction done per species pair, and the normalized score being fed into MCL program for the final gene family clustering; the other is different only by the way of obtaining the query-subject hits table via pair-wise diamond search in parallel instead of relying on the super cds file. The influential parameters are the same as in the whole paranome inference.

.. _collinearcoalescence:

Collinear coalescence inference of phylogeny
------------

Traditional phylogenetic inference relys on the sequence-based orthologues. The collinear orthologues are still understudied in inferring phylogeny. Such gene content- and order-conserved orthologues can be used to infer phylogeny under the multi-species coalescent (MSC) model with command below.

.. code-block:: console

   (ENV) $ wgd dmd cds1.fasta cds2.fasta cds3.fasta --anchorpoints apdata --segments smdata --listelements ledata --genetable gtdata --collinearcoalescence

We rely on the software ASTRAL to summarize the speceis tree based on the gene tree inferred per collinear orthologue group. At least 80% percent species is required to be present in the multiplicon and the intersected anchor points across all levels (the levels on the same scaffold are treated the same as on different scaffolds) are retreived as different collinear orthologue group (in the label for instance, Multiplicon1_Ap1, Multiplicon1_Ap2). The species occurance ratio can be controled by the parameter ``msogcut``. The gene tree inference method and parameter can be controlled by the parameter ``tree_method`` and ``treeset``.

.. py:function:: cli.dmd(sequences, outdir, tmpdir, cscore, inflation, eval, to_stop, cds, focus, anchorpoints, keepfasta, keepduplicates, globalmrbh, nthreads, orthoinfer, onlyortho, getnsog, tree_method, treeset, msogcut, geneassign, assign_method, seq2assign, fam2assign, concat, segments, listelements, collinearcoalescence, testsog, bins, buscosog, buscohmm, buscocutoff, genetable, normalizedpercent, nonormalization)

   Whole paranome inference

   Global/Local MRBHs inference

   Orthogroups inference

   Phylogeny inference based on the collinear coalescence

   :param sequences: Argument of sequence files.
   :type sequences: paths
   :param outdir: Path of desired output directory, default "wgd_dmd".
   :type outdir: str
   :param tmpdir: Path of temporary directory.
   :type tmpdir: str or None
   :param cscore: The c-score to restrict the homologs of MRBHs, default "None".
   :type cscore: float or None
   :param inflation: The inflation factor for MCL program, default "2.0".
   :type inflation: float
   :param eval: The e-value cut-off for similarity in diamond and/or hmmer, default "1e-10".
   :type eval: float
   :param to_stop: Whether to translate through STOP codons, default False.
   :type to_stop: boolean flags
   :param cds: Whether to only translate the complete CDS that starts with a valid start codon and only contains a single in frame stop codon at the end and must be dividable by three, default False.
   :type cds: boolean flags
   :param focus: The species to be merged on local MRBHs, default "None".
   :type focus: path or None
   :param anchorpoints: The anchor points data file, default "None".
   :type anchorpoints: path or None
   :param segments: The segments data file used in collinear coalescence analysis if initiated, default "None".
   :type segments: path or None
   :param listelements: The listsegments data file used in collinear coalescence analysis if initiated, default "None".
   :type listelements: path or None
   :param keepfasta: Whether to output the sequence information of MRBHs, default False.
   :type keepfasta: boolean flags
   :param keepduplicates: Whether to allow the same gene to occur in different MRBHs, default False.
   :type keepduplicates: boolean flags
   :param globalmrbh: Whether to initiate global MRBHs constructionhether to initiate global MRBHs construction, default False.
   :type globalmrbh: boolean flags
   :param nthreads: The number of threads to use, default 4.
   :type nthreads: int
   :param orthoinfer: Whether to initiate orthogroup infernece, default False.
   :type orthoinfer: boolean flags
   :param onlyortho: Whether to only conduct orthogroup infernece, default False.
   :type onlyortho: boolean flags
   :param getnsog: Whether to initiate the searching for nested single-copy gene families (NSOGs), default False.
   :type getnsog: boolean flags
   :param tree_method: which gene tree inference program to invoke, default "fasttree".
   :type tree_method: choice ['fasttree', 'iqtree', 'mrbayes']
   :param treeset: The parameters setting for gene tree inference, default None.
   :type treeset: multiple str options or None
   :param msogcut: The ratio cutoff for mostly single-copy family and species representation in collinear coalescence inference, default "0.8".
   :type msogcut: float
   :param geneassign: Whether to initiate the gene-to-family assignment analysis, default False.
   :type geneassign: boolean flags
   :param assign_method: Which method to conduct the gene-to-family assignment analysis, default "hmmer".
   :type assign_method: choice ['hmmer', 'diamond']
   :param seq2assign: The queried sequences data file in gene-to-family assignment analysis, default None.
   :type seq2assign: multiple path options or None
   :param fam2assign: The queried familiy data file in gene-to-family assignment analysis, default None.
   :type fam2assign: path or None
   :param concat: Whether to initiate the concatenation pipeline for orthogroup infernece, default False.
   :type concat: boolean flags
   :param collinearcoalescence: Whether to initiate the collinear coalescence analysis, default False.
   :type collinearcoalescence: boolean flags
   :param testsog: Whether to initiate the unbiased test of single-copy gene families, default False.
   :type testsog: boolean flags
   :param bins: The number of bins divided in gene length normalization, default "100".
   :type bins: int
   :param normalizedpercent: The percentage of upper hits used for gene length normalization, default "5".
   :type normalizedpercent: int
   :param nonormalization: Whether to call off the normalization, default False.
   :type nonormalization: boolean flags
   :param buscosog: Thether to initiate the busco-guided single-copy gene family analysis, default False.
   :type buscosog: boolean flags
   :param buscohmm: The HMM profile datafile in the busco-guided single-copy gene family analysis, default None.
   :type buscohmm: path or None
   :param buscocutoff: The HMM score cutoff datafile in the busco-guided single-copy gene family analysis, default None.
   :type buscocutoff: path or None
   :param genetable: The gene table datafile used in collinear coalescence analysis, default None.
   :type genetable: path or None

.. _ksdistribution:

*K*\ :sub:`S` distribution construction
------------

After obtaining the paralogous/orthologous gene family file, the construction of *K*\ :sub:`S` distribution can be achieved as below.

.. code-block:: console

   (ENV) $ wgd ksd families *cds.fasta

The gene family file is mandatory input, which can be acquired from the eariler steps using ``wgd dmd``. Depending on the number of species (or the provided cds sequence files), the meaning of constructed *K*\ :sub:`S` distribution differs. If only one species is given, the whole paranome *K*\ :sub:`S` distribution is to be established, which is a good material for the primary identification of WGDs. If orthogroups of multiple species are given, a paralog-ortholog mixed *K*\ :sub:`S` distribution is to be built that further subdivision per species pair of *K*\ :sub:`S` can be used in diagnosing the phylogenetic placement of focused WGD events. If the RBH gene family of two species is given, the constructed *K*\ :sub:`S` distribution is to show the *K*\ :sub:`S` age of divergence event, as illustrated in ``wgd v1``. This step is quite flexible in that various types of gene family files can be provided, for instance, the paralogous gene family, the orthogroup, the RBH gene family, the collinear orthologous gene family (inferred in the collinear coalescence inference analysis) and whatever gene families that users would like to calculate, as long as in the correct format (tab-separated table that each row represents a gene family while each column represents a species with the first column as the label of gene family and the first row as the label of species name). The influential parameters include the ``pairwise`` parameter determining whether to calculate the *K*\ :sub:`S` on the alignment basis of each paralogous gene pair instead of the whole alignment for that family, and the ``strip_gaps`` parameter controlling whether to drop all gaps in multiple sequence alignment, which only makes a difference when co-occurring with the ``pairwise`` parameter because the program ``codeml`` will strip all the gaps anyway. The default ``tree_method`` is the external software ``fasttree``, which can be replaced with the built-in ``cluster`` method.

.. _correctedksdistributio:

Corrected *K*\ :sub:`S` distribution construction
------------

To determine the phylogenetic location of a certain WGD, a rate-corrected *K*\ :sub:`S` distribution is needed to rescale the rate-dependent orthologous *K*\ :sub:`S` into the same rate of focused species. The command to achieve it is as below.

.. code-block:: console

   (ENV) $ wgd ksd families *cds.fasta --spair speciespair1 --spair speciespair2 --speciestree speciestree.nw


This step is mainly affected by the parameter ``reweight`` determining whether to recalculate the weight per species pair instead of using the weight calculated from the whole family, and the parameter ``onlyrootout`` determining whether to only conduct rate correction using the outgroup at root as outgroup. The parameter ``extraparanomeks`` enables the addition of extra paralogous *K*\ :sub:`S` data besides the probably existing ones in the original *K*\ :sub:`S` data and when duplicated paralogous gene pairs occur only the ones in extra paralogous *K*\ :sub:`S` data will be kept. This is to deal with the condition that users might provide gene families containing paralogous gene pairs which might not be complete and thus add more which might overlap with the extra paralogous *K*\ :sub:`S` data. But it's suggested to separate the paralogous *K*\ :sub:`S` only in the ``extraparanomeks`` parameter while the gene families only contain orthologous gene pairs (for instance, the global MRBHs) because the elmm modeling only considers the paralogs transmitted from the ``extraparanomeks`` parameter. Nonetheless, the ``extraparanomeks`` parameter would not affect the result of rate correction which only involves orthologous *K*\ :sub:`S`. We will discuss about it further in the recipe.

.. py:function:: cli.ksd(families, sequences, outdir, tmpdir, nthreads, to_stop, cds, pairwise, strip_gaps, tree_method,spair, speciestree, reweight, onlyrootout, extraparanomeks, anchorpoints, plotkde, plotapgmm, components)

   *K*\ :sub:`S` distribution construction

   Corrected *K*\ :sub:`S` distribution construction

   :param families: Argument of gene family file.
   :type families: path
   :param sequences: Argument of sequence files.
   :type sequences: paths
   :param outdir: Path of desired output directory, default "wgd_ksd".
   :type outdir: str
   :param tmpdir: Path of temporary directory.
   :type tmpdir: str or None
   :param nthreads: The number of threads to use, default 4.
   :type nthreads: int
   :param to_stop: Whether to translate through STOP codons, default False.
   :type to_stop: boolean flags
   :param cds: Whether to only translate the complete CDS that starts with a valid start codon and only contains a single in frame stop codon at the end and must be dividable by three, default False.
   :type cds: boolean flags
   :param pairwise: Whether to initiate pairwise *K*\ :sub:`S` estimation, default False.
   :type pairwise: boolean flag
   :param strip_gaps: Whether to drop all gaps in multiple sequence alignment.
   :type strip_gaps: boolean flag
   :param tree_method: which gene tree inference program to invoke, default "fasttree".
   :type tree_method: choice ['cluster', 'fasttree', 'iqtree']
   :param spair: The species pair to be plotted, default None.
   :type spair: multiple str options or None
   :param speciestree: The species tree to perform rate correction, default None.
   :type speciestree: path or None
   :param reweight: Whether to recalculate the weight per species pair, default False.
   :type reweight: boolean flag
   :param onlyrootout: Whether to only conduct rate correction using the outgroup at root as outgroup, default False.
   :type onlyrootout: boolean flag
   :param extraparanomeks: The extra paranome ks data to be plotted in the mixed *K*\ :sub:`S` distribution, default None.
   :type extraparanomeks: path or None
   :param anchorpoints: The anchorpoints.txt file to plot anchor *K*\ :sub:`S` in the mixed *K*\ :sub:`S` distribution, default None.
   :type anchorpoints: path or None
   :param plotkde: Whether to plot kde curve of orthologous *K*\ :sub:`S` distribution over histogram in the mixed Ks distribution, default False.
   :type plotkde: boolean flag
   :param plotapgmm: Whether to plot mixture modeling of anchor *K*\ :sub:`S` in the mixed *K*\ :sub:`S` distribution, default False.
   :type plotapgmm: boolean flag
   :param plotelmm: Whether to plot elmm mixture modeling of paranome *K*\ :sub:`S` in the mixed *K*\ :sub:`S` distribution, default False.
   :type plotelmm: boolean flag
   :param components: The range of the number of components to fit in anchor *K*\ :sub:`S` mixture modeling, default (1,4).
   :type components: (int,int)

.. _mixture:

Mixture model clustering of *K*\ :sub:`S` distribution
------------

This function inherits from ``wgd v1`` that ``bgmm`` and ``gmm`` are available for mixture modeling upon any *K*\ :sub:`S` distribution. The command to achieve it is as below.

.. code-block:: console

   (ENV) $ wgd mix ks.tsv

This part of analysis is mainly relying on the mixture module of scikit-learn library. Almost all the parameters of this function will have an impact on the results. Please see the description below for the detailed information.

.. py:function:: cli.mix(ks_distribution, filters, ks_range, method, components, bins, output_dir, gamma, n_init, max_iter)

   Mixture modeling of *K*\ :sub:`S` distribution

   :param ks_distribution: Argument of *K*\ :sub:`S` distribution file.
   :type ks_distribution: path
   :param filters: The cutoff alignment length, default "300".
   :type filters: int
   :param ks_range: The *K*\ :sub:`S` range to be considered, default (0, 5)
   :type ks_range: (float,float)
   :param method: which mixture model to use, default "gmm".
   :type method: choice ['gmm', 'bgmm']
   :param components: The range of the number of components to fit, default (1,4)
   :type components: (int,int)
   :param bins: The number of bins in *K*\ :sub:`S` distribution, default "50".
   :type bins: int
   :param outdir: Path of desired output directory, default "wgd_mix".
   :type outdir: str
   :param gamma: The gamma parameter for bgmm models, default "0.001".
   :type gamma: float
   :param n_init: The number of k-means initializations, default "200".
   :type n_init: int
   :param max_iter: The maximum number of iterations, dafault "1000".
   :type max_iter: int

.. _synteny:

Synteny inference
------------

Synteny is frequently called nowadays in profiling the evolution of psedochromosomes. The program ``wgd syn`` does all the inference work about synteny. The simplest command is as below.

.. code-block:: console

   (ENV) $ wgd syn families.tsv gff3

The influential parameters for synteny inference include the ``minlen`` controling the minimum length of a scaffold to be considered, the ``maxsize`` controling the maximum family size to be considered, the ``minseglen`` determining the minimum length of segments to considered and the ``keepredun`` determining whether to keep redundant multiplicons.

.. py:function:: cli.syn(families, gff_files, ks_distribution, outdir, feature, attribute, minlen, maxsize, ks_range, iadhore_options, ancestor, minseglen, keepredun)

   Synteny inference

   :param families: Argument of gene family file.
   :type families: path
   :param gff_files: Argument of gff3 files.
   :type gff_files: paths
   :param ks_distribution: The *K*\ :sub:`S` distribution datafile, default None.
   :type ks_distribution: path or None
   :param outdir: Path of desired output directory, default "wgd_syn".
   :type outdir: str
   :param feature: The feature for parsing gene IDs from GFF files, default "gene".
   :type feature: str
   :param attribute: The attribute for parsing the gene IDs from the GFF files, default "ID".
   :type attribute: str
   :param minlen: The minimum length of a scaffold to be included in dotplot, default "-1".
   :type minlen: int
   :param maxsize: The maximum family size to include, default "200".
   :type maxsize: int
   :param ks_range: The *K*\ :sub:`S` range in colored dotplot, default (0, 5).
   :type ks_range: (float,float)
   :param iadhore_options: The parameter setting in iadhore, default "".
   :type iadhore_options: str
   :param ancestor: The assumed ancestor species, it's still under development, default None.
   :type ancestor: str or None
   :param minseglen: The minimum length of segments to include, ratio if <= 1, default 100000.
   :type minseglen: float
   :param keepredun: Whether to keep redundant multiplicons, default False.
   :type keepredun: boolean flags

.. _searchanchor:

Search of anchor pairs for molecular dating
------------

If we want to date the identified WGD events, determining anchor pairs to be dated is the key step. The program ``wgd peak`` can achieve this goal. The command is as below.

.. code-block:: console

   (ENV) $ wgd peak ksdata -ap apdata -sm smdata -le ledata -mp mpdata

There are two methods which can be called in this step. One is the heuristic method which can be called by adding the flag ``heuristic``, where a heuristic refinement of anchor pairs based on the detected peaks from ``scipy.signal`` and the associated properties is conducted. Another method is to retreive the highest density region (HDR) of the segment-guided syntelogs with users-defined *K*\ :sub:`S` saturation cutoff.

.. py:function:: cli.peak(ks_distribution, anchorpoints, outdir, alignfilter, ksrange, bin_width, weights_outliers_included, method, seed, em_iter, n_init, components, boots, weighted, plot, bw_method, n_medoids, kdemethod, n_clusters, kmedoids, guide, prominence_cutoff, kstodate, family, rel_height, ci, manualset, segments, hdr, heuristic, listelements, multipliconpairs, kscutoff, gamma)

   Search of anchor pairs used in WGD dating

   :param ks_distribution: Argument of *K*\ :sub:`S` distribution datafile.
   :type ks_distribution: path
   :param anchorpoints: The anchor points datafile, default None.
   :type anchorpoints: path or None
   :param segments: The segments datafile, default None.
   :type segments: path or None
   :param outdir: Path of desired output directory, default "wgd_peak".
   :type outdir: str
   :param alignfilter: The cutoff for alignment identity, length and coverage, default (0.0, 0, 0.0).
   :type alignfilter: (float,int,float)
   :param ksrange: The range of *K*\ :sub:`S` to be analyzed, default (0,5).
   :type ksrange: (float,float)
   :param bin_width: The bandwidth of *K*\ :sub:`S` distribution, default "0.1".
   :type bin_width: float
   :param weights_outliers_included: Whether to include *K*\ :sub:`S` outliers, default False.
   :type weights_outliers_included: boolean flags
   :param method: Which mixture model to use, default "gmm".
   :type method: choice ['gmm', 'bgmm']
   :param seed: Random seed given to initialization, default "2352890".
   :type seed: int
   :param em_iter: The number of EM iterations to perform, default "200".
   :type em_iter: int
   :param n_init: The number of k-means initializations, default "200".
   :type n_init: int
   :param components: The range of the number of components to fit, default (1,4).
   :type components: (int,int)
   :param boots: The number of bootstrap replicates of kde, default "200".
   :type boots: int
   :param weighted: Whether to use node-weighted method of de-redundancy, default False.
   :type weighted: boolean flags
   :param plot: The plotting method to be used, default "identical".
   :type plot: choice(['stacked', 'identical'])
   :param bw_method: The bandwidth method to be used, default "silverman".
   :type bw_method: choice['silverman', 'ISJ']
   :param n_medoids: The number of medoids to fit, default "2".
   :type n_medoids: int
   :param kdemethod: The kde method to be used, default "scipy".
   :type kdemethod: choice['scipy', 'naivekde', 'treekde', 'fftkde']
   :param n_clusters: The number of clusters to plot Elbow loss function, default "5".
   :type n_clusters: int
   :param kmedoids: Whether to initiate K-Medoids clustering analysis, default False.
   :type kmedoids: boolean flags
   :param guide: The regime residing anchors, default "segment".
   :type guide: choice['multiplicon', 'basecluster', 'segment']
   :param prominence_cutoff: The prominence cutoff of acceptable peaks, default "0.1".
   :type prominence_cutoff: float
   :param kstodate: The range of *K*\ :sub:`S` to be dated in heuristic search, default (0.5, 1.5).
   :type kstodate: (float,float)
   :param family: The family to filter *K*\ :sub:`S` upon, default None.
   :type family: path or None
   :param manualset: Whether to output anchor pairs with manually set *K*\ :sub:`S` range, default False.
   :type manualset: boolean flags
   :param rel_height: The relative height at which the peak width is measured, default "0.4".
   :type rel_height: float
   :param ci: The confidence level of log-normal distribution to date, default "95".
   :type ci: int
   :param hdr: The highest density region (HDR) in a given distribution to date, default "95".
   :type hdr: int
   :param heuristic: Whether to initiate heuristic method of defining CI for dating, default False.
   :type heuristic: boolean flags
   :param listelements: The listsegments datafile, default None.
   :type listelements: path or None
   :param multipliconpairs: The multipliconpairs datafile, default None.
   :type multipliconpairs: path or None
   :param kscutoff: The *K*\ :sub:`S` saturation cutoff in dating, default "5".
   :type kscutoff: float
   :param gamma: The gamma parameter for bgmm models, default "1e-3".
   :type gamma: float

.. _phylogeneticinference:

Concatenation-/Coalescence-based phylogenetic inference
------------

The program ``wgd focus`` can handle various analysis. To recover the phylogeny, either in concatenation- or coalescence-based method, the following command can be used.

.. code-block:: console

   (ENV) $ wgd focus families *cds.fasta --concatenation/--coalescence

What this step does is 1) write sequence for each gene family, 2) infer multiple sequence alignment (MSA) for each gene family, if with concatenation method, 3) concatenate all MSA and infer the maximum likelihood gene tree as species tree, if with coalescence method 3) infer maximum likelihood gene tree for each gene family and summarize the species tree using ASTRAL. The influential parameters include the ``tree_method`` parameter controling the program to infer gene tree and the ``treeset`` parameter to control the parameter setting for gene tree inference.

.. _annotation:

Functional annotation of gene families
------------

The functional annotation of gene families can also be achieved in ``wgd focue`` using the command below.

.. code-block:: console

   (ENV) $ wgd focus families *cds.fasta --annotation eggnog --eggnogdata eddata --dmnb dbdata

This step relies on the database provided by users. It's required to pre-set the executable of annotation program for instance, ``eggnog``, ``hmmer`` and ``interproscan`` in the environment variable. The influential parameter is the ``evalue`` that controls the e-value cut-off for annotation.

.. _dating:

WGD dating
------------

The program ``wgd focus`` is the final step in WGD dating. After obtaining the anchor pairs-merged orthogroups from ``wgd dmd``, the phylogenetic dating can be conducted with the command below.

.. code-block:: console

   (ENV) $ wgd focus families sequence1 sequence2 sequence3 --dating mcmctree --speciestree species_tree.nw

The command shown above is a simple example that calls the molecular dating program ``mcmctree`` to conduct the dating. A species tree is mandatory input. Further dating parameters can be provided with ``datingset`` parameters.

.. py:function:: cli.focus(families, sequences, outdir, tmpdir, nthreads, to_stop, cds, strip_gaps, aligner, tree_method, treeset, concatenation, coalescence, speciestree, dating, datingset, nsites, outgroup, partition, aamodel, ks, annotation, pairwise, eggnogdata, pfam, dmnb, hmm, evalue, exepath, fossil, rootheight, chainset, beastlgjar, beagle, protdating)

   Concatenation-/Coalescence-based phylogenetic inference
   Functional annotation of gene families
   Phylogenetic dating of WGDs

   :param families: Argument of gene family file.
   :type families: path
   :param sequences: Argument of sequence files.
   :type sequences: paths
   :param outdir: Path of desired output directory, default "wgd_focus".
   :type outdir: str
   :param tmpdir: Path of temporary directory.
   :type tmpdir: str or None
   :param nthreads: The number of threads to use, default 4.
   :type nthreads: int
   :param to_stop: Whether to translate through STOP codons, default False.
   :type to_stop: boolean flags
   :param cds: Whether to only translate the complete CDS that starts with a valid start codon and only contains a single in frame stop codon at the end and must be dividable by three, default False.
   :type cds: boolean flags
   :param strip_gaps: Whether to drop all gaps in multiple sequence alignment.
   :type strip_gaps: boolean flag
   :param aligner: Which alignment program to use, default "mafft".
   :type aligner: choice(['muscle', 'prank', 'mafft'])
   :param tree_method: which gene tree inference program to invoke, default "fasttree".
   :type tree_method: choice ['fasttree', 'iqtree', 'mrbayes']
   :param treeset: The parameters setting for gene tree inference, default None.
   :type treeset: multiple str options or None
   :param concatenation: Whether to initiate the concatenation-based species tree inference, default False.
   :type concatenation: boolean flag
   :param coalescence: Whether to initiate the coalescence-based species tree inference, default False.
   :type coalescence: boolean flag
   :param speciestree: The species tree datafile for dating, default None.
   :type speciestree: path or None
   :param dating: Which molecular dating program to use, default None.
   :type dating: choice(['beast', 'mcmctree', 'r8s', 'none'])
   :param datingset: The parameters setting for dating program, default None.
   :type datingset: multiple str options or None
   :param nsites: The nsites information for r8s dating, default None.
   :type nsites: int
   :param outgroup: The outgroup information for r8s dating, default None.
   :type outgroup: str
   :param partition: Whether to initiate partition dating analysis for codon, default False.
   :type partition: boolean flag
   :param aamodel: Which protein model to be used in mcmctree, default "poisson".
   :type aamodel: choice(['poisson','wag', 'lg', 'dayhoff'])
   :param ks: Whether to initiate *K*\ :sub:`S` calculation, default False.
   :type ks: boolean flag
   :param annotation: Which annotation program to use, default None.
   :type annotation: choice(['none','eggnog', 'hmmpfam', 'interproscan'])
   :param pairwise: Whether to initiate pairwise *K*\ :sub:`S` estimation, default False.
   :type pairwise: boolean flag
   :param eggnogdata: The eggnog annotation datafile, default None.
   :type eggnogdata: path or None
   :param pfam: Which option to use for pfam annotation, default None.
   :type pfam: choice(['none', 'denovo', 'realign'])
   :param dmnb: The diamond database for annotation, default None.
   :type dmnb: path or None
   :param hmm: The HMM profile for annotation, default None.
   :type hmm: path or None
   :param evalue: The e-value cut-off for annotation, default "1e-10".
   :type evalue: float
   :param exepath: The path to the interproscan executable, default None.
   :type exepath: path or None
   :param fossil: The fossil calibration information in Beast, default ('clade1;clade2', 'taxa1,taxa2;taxa3,taxa4', '4;5', '0.5;0.6', '400;500').
   :type fossil: (str,str,str,str,str)
   :param rootheight: The root height calibration info in Beast, default (4,0.5,400).
   :type rootheight: (float,float,float)
   :param chainset: The parameters of MCMC chain in Beast, default (10000,100).
   :type chainset: (int,int)
   :param beastlgjar: The path to beastLG.jar, default None.
   :type beastlgjar: path or None
   :param beagle: Whether to use beagle in Beast, default False.
   :type beagle: boolean flag
   :param protdating: Whether to only initiate the protein-concatenation based dating analysis, default False.
   :type protdating: boolean flag

.. _viz:

*K*\ :sub:`S` distribution visualization
------------

The program ``wgd viz`` can be used in plotting *K*\ :sub:`S` distribution and synteny. To visualize the *K*\ :sub:`S` distribution, the command below can be used.

.. code-block:: console

   (ENV) $ wgd viz --data ks.tsv

The program ``wgd viz`` will automately calculate and plot the exponential-lognormal mixture modeling (ELMM) result. The influential parameters include the ``em_iterations`` controling the maximum EM iterations and the ``em_initializations`` controling the the maximum EM initializations, the ``prominence_cutoff`` determining the prominence cutoff of acceptable peaks and the ``rel_height`` to set the relative height at which the peak width is measured.

Corrected *K*\ :sub:`S` distribution visualization
------------

To add rate correction result into the *K*\ :sub:`S` plot, one can use the command below.

.. code-block:: console

   (ENV) $ wgd viz --data ks.tsv --spair speciespair1 --spair speciespair2 --speciestree speciestree.nw --gsmap gene_species.map

The additional required file ``gene_species.map`` is automately produced from the ``wgd ksd`` step when producing the ``ks.tsv`` file.

Synteny visualization
------------

The synteny plot produced by the program ``wgd syn`` can be reproduced by ``wgd viz`` too. The command is as below.

.. code-block:: console

   (ENV) $ wgd viz --anchorpoints apdata --segments smdata --multiplicon mtdata --genetable gtdata

The influential parameters include the ``minlen`` controling the minimum length of a scaffold to be included in dotplot, the ``maxsize`` determining the maximum family size to include, the ``minseglen`` determining the minimum length of segments to include, the ``keepredun`` controling whether to keep redundant multiplicons.

.. py:function:: cli.viz(datafile,spair,outdir,gsmap,plotkde,reweight,em_iterations,em_initializations,prominence_cutoff,segments,minlen,maxsize,anchorpoints,multiplicon,genetable,rel_height,speciestree,onlyrootout,minseglen,keepredun,extraparanomeks,plotapgmm,plotelmm,components)

   *K*\ :sub:`S` distribution visualization
   Synteny visualization

   :param datafile: The path to datafile, default None.
   :type datafile: path or None
   :param spair: The species pair to be plotted, default None.
   :type spair: multiple str options
   :param outdir: Path of desired output directory, default "wgd_focus".
   :type outdir: str
   :param gsmap: The gene name-species name map, default None.
   :type gsmap: path or None
   :param plotkde: Whether to plot kde curve upon histogram, default False.
   :type plotkde: boolean flag
   :param reweight: Whether to recalculate the weight per species pair, default False.
   :type reweight: boolean flag
   :param em_iterations: The maximum EM iterations, default "200".
   :type em_iterations: int
   :param em_initializations: The maximum EM initializations, default "200".
   :type em_initializations: int
   :param prominence_cutoff: The prominence cutoff of acceptable peaks, default "0.1".
   :type prominence_cutoff: float
   :param segments: The segments data file, default None.
   :type segments: path or None
   :param minlen: The minimum length of a scaffold to be included in dotplot, if "-1" was set, the 10% of the longest scaffolds will be used, default "-1".
   :type minlen: int
   :param maxsize: The maximum family size to include, default "200".
   :type maxsize: int
   :param anchorpoints: The anchor points datafile, default None.
   :type anchorpoints: path or None
   :param multiplicon: The multiplicons datafile, default None.
   :type multiplicon: path or None
   :param genetable: The gene table datafile, default None.
   :type genetable: path or None
   :param rel_height: The relative height at which the peak width is measured, default "0.4".
   :type rel_height: float
   :param speciestree: The species tree to perform rate correction, default None.
   :type speciestree: path or None
   :param onlyrootout: Whether to only conduct rate correction using the outgroup at root as outgroup, default False.
   :type onlyrootout: boolean flag
   :param minseglen: The minimum length of segments to include in ratio if <= 1, default "100000".
   :type minseglen: float
   :param keepredun: Whether to keep redundant multiplicons, default False.
   :type keepredun: boolean flag
   :param extraparanomeks: The extra paranome ks data to be plotted in the mixed *K*\ :sub:`S` distribution, default None.
   :type extraparanomeks: path or None
   :param plotapgmm: Whether to plot mixture modeling of anchor *K*\ :sub:`S` in the mixed *K*\ :sub:`S` distribution, default False.
   :type plotapgmm: boolean flag
   :param plotelmm: Whether to plot elmm mixture modeling of paranome *K*\ :sub:`S` in the mixed *K*\ :sub:`S` distribution, default False.
   :type plotelmm: boolean flag
   :param components: The range of the number of components to fit in anchor *K*\ :sub:`S` mixture modeling, default (1,4).
   :type components: (int,int) 

