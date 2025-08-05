#import "@preview/peace-of-posters:0.5.6" as pop
#import "@preview/xarrow:0.3.0": xarrow

#let theme = (
    "body-box-args": (
        inset: (x: 0.0em, y: 0.6em),
        width: 100%,
        stroke: none,
    ),
    "body-text-args": (:),
    "heading-box-args": (
        inset: 0em,
        width: 100%,
        stroke: none,
    ),
    "heading-text-args": (
        fill: rgb("#008b92"),
        weight: "bold",
    ),
)

#set page("a0", margin: 3cm)
#pop.set-poster-layout(pop.layout-a0)
#pop.set-theme(theme)
#set text(font: "Futura")
#let box-spacing = 1.5em
#set columns(gutter: box-spacing)
#set block(spacing: box-spacing)
#pop.update-poster-layout(
  spacing: box-spacing,
  heading-size : 40pt,
  title-size	: 80pt,
  //body-size: 33pt,
  authors-size: 47pt,
  institutes-size: 33pt,
)

#pop.title-box(
    [
        #set text(fill: white)
        //#image("eveco-en-white.svg", width: 25%)
        Model-based demographic inference of recent invasions from genomic data
    ],
    authors: [
        //#v(0.5cm)
        #set text(fill: black)
        Francisco Campuzano Jiménez#super("1"),
        Els De Keyzer#super("1"),
        Arthur Zwaenepoel#super("1"),
        Hannes Svardal#super("1,2")
    ],
    institutes: [
        #set text(fill: black, weight: "regular")
        #super("1")Evolutionary Ecology Group, University of Antwerp, Belgium
        #super("2")Naturalis Biodiversity Center, Leiden, Netherlands
    ],
  //background: box(image("pink-yellow.png", height: 16cm, width: 100%), inset: -2cm, outset: 0pt),
  background: box(rect(fill: rgb("#008b92"), width: 100%, height: 13cm), inset: -3cm, outset: 0pt),
)
#v(1cm)
#grid(
    gutter: pop.layout-a0.spacing * 4,
    column-gutter: pop.layout-a0.spacing * 4,
    columns: (55%, 45%),
    box(width: 100%)[
      #pop.column-box(heading: "Introduction")[
        #v(0.5cm)
        - Biological invasions are a major threat to biodiversity but also can be viewed as large-scale, unplanned experiments on how populations adapt to novel environments.
        - Accurate estimates of recent effective population size ($N_e$) can guide conservation management plans and help disentangle the effects of selection and demography on genetic variation.
        - Existing methods to infer very recent effective population size are not well-suited for realistic invasion scenarios.
      ],
      #pop.column-box(heading: "Theoretical framework")[
        #v(0.5cm)
        We propose a method to infer recent evolution in $N_e$ that can incorporate relevant prior information of the invasion process and the source population.
        #image("combined_figure_example.svg", width: 100%)
      ]
      #pop.column-box(heading: "Modelling multiple introductions")[
      #v(0.5cm)
      Multiple introductions are common in biological invasions, although difficult to model. In order to accommodate them, we leveraged recent advances in inhomogeneous phase-type distributions to compute approximate likelihoods.
         #figure(
           image("multiple_introductions.svg", width: 100%),
           caption: text(
             [_Figure 1: Likelihoods were obtained via numerical integration of the cumulative distribution function of the TMCRA using `phasegen`, keeping all other parameters fixed at their true values._],
             size: 26pt,
           ), numbering: none
         )
         Preliminary results suggest we can estimate the migration rate $m$ from LD data of the focal population solely. However, the posterior distribution for all other parameters is highly multimodal and raises concerns about identifiability.
      ],
    ],
    box(width: 100%)[
    #pop.column-box(heading: "Method comparison")[
    #v(0.5cm)
      The proposed method works better than existing state-of-the-art methods, at least in part because it does not penalize an abrupt change in $N_e$ at the beginning of the invasion.
    #figure(
      image("method_boxplot.svg", width: 90%),
      caption: text(
        [_Figure 2: We computed the absolute relative error of $N_e$ across the last 100 generations in 9 demographic scenarios (25 replicates each). Each dataset consisted of 200 diploid individuals and 25 chromosomes. HapNe-IBD has higher error than other methods (not shown)._ ],
        size: 26pt
      ),
      numbering: none
    )]
    #pop.column-box(heading: "Estimating key parameters")[
    #v(0.5cm)
    The time of invasion and the order of magnitude of the effective founder size can be estimated precisely, even without any _a priori_ knowledge of the source population.
    #figure(
      image("expected_posterior.svg", width: 100%),
      caption: text(
        [_Figure 3: We compared the expected posterior with the ground truth of the effective founder size and time of invasion across 225 synthetic datasets. Doing inference, the prior of $N_e(t)$ was chosen to match observed genetic diversity._],
        size: 26pt
      ),
      numbering: none
    )
    ],
    #pop.column-box(heading: "Take home message")[
      #v(0.5cm)
      - Model-based methods based on LD data can provide better estimates of demographic parameters for invasion scenarios than existing methods.
      - Recent advances in computational and statistical methods enable modeling more complex demographic scenarios and going beyond overly simplistic models.
      ],
      #v(1cm)
      #figure(
        image("qrcode.svg", width: 50%),
        numbering: none
      )
      #set align(center)
      #image("EVECO_customlogo.png", width: 70%)
    ]


)


#v(-1cm)
#pop.bottom-box(
    heading-box-args: (inset: 0.5cm, fill: rgb("#008b92")),
)[
  #box(height: 2em)[#set text(size: 33pt, fill: white, )
    //#align(horizon)[#text("@currocam")  #h(1fr) #text("currocam.github.io")]]
    #align(horizon)[#text("You can find me on social media as @currocam") #h(1fr) #text("currocam.github.io") #h(1fr) #text("curro.campuzanojimenez@uantwerpen.be")]]

]
