---
layout: post
title: "Thermo- and Mechanical Property Prediction of Amorphous Cellulose Triacetate Plastic"
author: "Samuel Holmes"
---
Hello, and welcome to my blog! I’m Sam, a computational chemist and postdoctoral researcher in the Woods Lab at the CCRC at the University of Georgia, and I specialize in force field parameter development and molecular modeling of glycomaterials. This blog documents my ongoing efforts to model and design eco-friendly polysaccharide-based plastics using all-atom MD.

Thermoplastics have become some of the most widely used and indispensable materials in the modern world. They are used in everything from medical implants and drug delivery systems to toys, textiles, packaging, automotive parts, building materials, and even military-grade composites. As global plastic consumption continues to rise with a growing population, the need for renewable, sustainable alternatives to petroleum-based plastics has become more urgent than ever. 

Fortunately, natural and engineered polysaccharides offer a vast and versatile feedstock for designing next-generation thermoplastics. Some have already demonstrated real-world potential, like cellulose triacetate (CTA), which was used in the production of LEGO bricks for decades. However, synthesizing novel glycomaterials with defined chemical structures remains a significant challenge. This is mainly due to the high density of chemically similar hydroxyl groups along the polysaccharide backbone, which complicates selective substitution using traditional protecting group strategies. These synthetic difficulties have motivated the development of computational approaches that can predict material properties as a function of chemical structure, helping to narrow the design space and accelerate the discovery of sustainable glycomaterial thermoplastics.

To that end, in this blog series, I use all-atom MD simulations with the GLYCAM06 force field in AMBER to predict thermo- and mechanical properties of an amorphous cellulose triacetate plastic. This is a living research project and will evolve over several months—expect drafts, breakthroughs, and the occasional setback. Thanks for joining me.


# Goals
- Successfully build and parameterize a system of 20 x CTA with DP=100.
- Equilibrate the system at 600 K to allow chain entanglement, and then cool down to 300 K to achieve a "realistic" starting configuration of an amorphous CTA plastic.
- Predict glass transition temperature and Young's modulus. 

## Part 1: Building and Parameterizing the System

One of the first decisions I had to make in modeling cellulose triacetate (CTA) was how closely to match the degree of substitution (DS) in real CTA plastics. The DS refers to the average number of hydroxyl groups per glucose unit that have been replaced with acetate groups, with a maximum of 3.0. In practice, polymer chemists rarely achieve a perfect DS of exactly 3.0. Even when aiming for full substitution, the final product typically falls somewhere in the range of 2.7 to just under 3.0, depending on reaction conditions and steric accessibility of the hydroxyls.

These small differences in DS can affect solubility, crystallinity, and mechanical properties, so it’s important to at least think critically about how we represent them in simulations. However, modeling a statistically accurate distribution of DS across a long polymer chain would introduce a lot of complexity—especially for an all-atom force field like GLYCAM that relies on clearly defined residue types and parameterized chemical environments.

So for this initial model, I decided to use a simplified approach: a fully substituted CTA with DS = 3.0. Every glucose unit is triacetylated. While not perfectly reflective of a real-world CTA sample, this version is chemically well-defined, force field-compatible, and likely close enough in structure to exhibit similar bulk behavior—especially in the amorphous phase. If needed, I can revisit partial substitution in future models using a randomized or pattern-based DS distribution.
