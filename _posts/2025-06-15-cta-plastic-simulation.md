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

## Part 1: Building and Parameterizing the Initial System

The first objective in this project is to obtain a 3D model of CTA that is physically realistic and consistent with the DS of real CTA plastics. The DS refers to the average number of hydroxyl groups per glucose unit that have been replaced with acetate groups, with a maximum of 3.0. In practice, polymer chemists rarely achieve a perfect DS of exactly 3.0. Even when aiming for full substitution, the final product typically falls somewhere in the range of 2.7 to just under 3.0, depending on reaction conditions and steric accessibility of the hydroxyls. These small differences in DS can affect the the thermoplastic properties of CTA, which is why it is important to try and match the real DS. 

For this initial model, I am going to use a simplified approach: a fully substituted CTA with DS = 3.0, meaning that every glucose monomer is acetyated at the 2-, 3-, and 6- positions. While not perfectly reflective of a real-world CTA sample, this version is chemically well-defined, and likely close enough in chemical structure to exhibit similar bulk behavior—especially in the amorphous phase. If needed, I can revisit partial substitution in future models using a randomized or pattern-based DS distribution.

### Obtaining a 3D Model of a Single CTA 100-mer (DS = 3.0)

In this project, I will use the Glycam Molecular Modeling Library (GMML) to construct a 3D model of CTA, which is a C++ software package written by Oliver Grant, a research scientist in the Woods Lab. I've included Oliver's GitHub URL here if you are interested in checking out his coding projects (https://github.com/gitoliver). A new version of GMML has recently been released to the public, called GMML2. GMML2 is an incredibly valuable piece of software, as it was prepared with the strict formatting requirements of LEaP in mind - the notoriously finnicky program in AMBER that attempts to match force field parameters to a 3D model. In short, GMML2 is an ideal choice here for model generation as it prepares a physically realistic model in PDB format with all the necessary information for LEaP to correctly parameterize it, in particular, GLYCAM atom and residue names, TER cards to separate each residue, and bonding information. GMML2 can be installed as a standalone package and used on the command line, and information about installation and usage is included here (https://github.com/GLYCAM-Web/gmml2). 

To build the model for CTA, we can take advantage of some helpful GMML functionality called Glycam-Condensed Notation, which allows glycomaterials to be built from simple strings. More information can be found here (https://glycam.org/docs/custombuilders/condensed-notation/index.html). 

Once GMM2 has be succesfully installed, we can build the model of CTA with the carbohydrate builder:

carbohydrateBuilder [input file] _ [output directory]

This command takes three arguments:
1) You input file containing the glycomaterial you want to build.
2) The delimiter used in your input file.
3) The output directory for your model.

3) is straightforward, but let's break down 1) and 2). My input file looks like this:

CTA.100mer_[4DGlcp[2Ac,3Ac,6Ac]b1-]<100>OH

The format the carbohydrate builder is expecting is [sequence_name][delimiter][sequence in glycam condensed notation]. You can choose whatever delimeter you want, but the second argument in the command must specify what the delimeter is.




