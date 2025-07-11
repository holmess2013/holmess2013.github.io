---
layout: post
title: "Thermo- and Mechanical Property Prediction of Amorphous Cellulose Triacetate Plastic"
author: "Samuel Holmes"
tags: [MD-simulation, cellulose, force-field, thermoplastic, AMBER, GLYCAM, glycomaterials]

---
# Introduction
Hello, and welcome to my blog! I’m Sam, a computational chemist and postdoctoral researcher in the Woods Lab at the CCRC at the University of Georgia, and I specialize in force field parameter development and molecular modeling of glycomaterials. This blog documents my ongoing efforts to model and design eco-friendly polysaccharide-based plastics using all-atom MD.

Thermoplastics have become some of the most widely used and indispensable materials in the modern world. They are used in everything from medical implants and drug delivery systems to toys, textiles, packaging, automotive parts, building materials, and even military-grade composites. As global plastic consumption continues to rise with a growing population, the need for renewable, sustainable alternatives to petroleum-based plastics has become more urgent than ever. 

Fortunately, natural and engineered polysaccharides offer a vast and versatile feedstock for designing next-generation thermoplastics. Some have already demonstrated real-world potential, like cellulose triacetate (CTA), which was used in the production of LEGO bricks for decades. However, synthesizing novel glycomaterials with defined chemical structures remains a significant challenge. This is mainly due to the high density of chemically similar hydroxyl groups along the polysaccharide backbone, which complicates selective substitution using traditional protecting group strategies. These synthetic difficulties have motivated the development of computational approaches that can predict material properties as a function of chemical structure, helping to narrow the design space and accelerate the discovery of sustainable glycomaterial thermoplastics.

To that end, in this blog series, I use all-atom MD simulations with the GLYCAM06 force field in AMBER to predict thermo- and mechanical properties of an amorphous cellulose triacetate plastic. This is a living research project and will evolve over several months—expect drafts, breakthroughs, and the occasional setback. Thanks for joining me.


# Goals
- Successfully build and parameterize a system of 20 x CTA polymers with DP = 100 and DS = 3.0.
- Equilibrate the system at 600 K to allow chain entanglement, and then cool down to 300 K to achieve a "realistic" starting configuration of an amorphous CTA plastic.
- Predict glass transition temperature and Young's modulus. 

## Part 1: Building and Parameterizing the Initial System

The first objective in this project is to obtain a 3D model of CTA that is physically realistic and consistent with the DS of real CTA plastics. The DS refers to the average number of hydroxyl groups per glucose unit that have been replaced with acetate groups, with a maximum of 3.0. In practice, polymer chemists rarely achieve a perfect DS of exactly 3.0. Even when aiming for full substitution, the final product typically falls somewhere in the range of 2.7 to just under 3.0, depending on reaction conditions and steric accessibility of the hydroxyls. These small differences in DS can affect the the thermoplastic properties of CTA, which is why it is important to try and match the real DS. 

For this initial model, I am going to use a simplified approach: a fully substituted CTA with DS = 3.0, meaning that every glucose monomer is acetyated at the 2-, 3-, and 6- positions. While not perfectly reflective of a real-world CTA sample, this version is chemically well-defined, and likely close enough in chemical structure to exhibit similar bulk behavior—especially in the amorphous phase. If needed, I can revisit partial substitution in future models using a randomized or pattern-based DS distribution.

### Obtaining a 3D Model of a Single CTA 100-mer (DS = 3.0)

In this project, I will use the Glycam Molecular Modeling Library (GMML) to construct a 3D model of CTA, which is a C++ software package written by Oliver Grant, a research scientist in the Woods Lab. I've included Oliver's GitHub URL here if you are interested in checking out his coding projects (https://github.com/gitoliver). A new version of GMML has recently been released to the public, called GMML2. GMML2 is an incredibly valuable piece of software, as it was prepared with the strict formatting requirements of LEaP in mind - the notoriously finnicky program in AMBER that attempts to match force field parameters to a 3D model. In short, GMML2 is an ideal choice here for model generation as it rapidly generates physically realistic molecular models of glycomaterials in PDB format with all the necessary information for LEaP to correctly parameterize it, in particular, GLYCAM atom and residue names, TER cards to separate each residue, and bonding information. GMML2 can be installed as a standalone package and used on the command line, and information about installation and usage is included here (https://github.com/GLYCAM-Web/gmml2). 

To build the model for CTA, we can take advantage of some helpful GMML functionality called Glycam-Condensed Notation, which allows glycomaterials to be built from simple strings. More information can be found here (https://glycam.org/docs/custombuilders/condensed-notation/index.html). 

Once GMML2 has be succesfully installed, we can build the model of CTA with the carbohydrate builder:

```
carbohydrateBuilder [input file] _ [output directory]
```

This command takes three arguments:
1) You input file containing the glycomaterial you want to build.
2) The delimiter used in your input file.
3) The output directory for your model.

3) is straightforward, but let's break down 1) and 2). My input file looks like this:
```
CTA.100mer_[4DGlcp[2Ac,3Ac,6Ac]b1-]<100>OH
```
The format the carbohydrate builder is expecting is [sequence_name][delimiter][sequence in glycam condensed notation]. You can choose whatever delimeter you want, but the second argument in the command must specify what the delimeter is.

After building the sequence, GMML2 generates two files, the pdb containing your structure and an AMBER off file. We can ignore the latter. Let's visualize the results of the build in VMD:

![GMML2_generated_CTA_model](/figures/GMML_CTA_figure1.png)

As you can see, its a pretty massive, linear polysaccharide, so I show both the entire 100-mer as well as a zoomed-in view so that you can see the molecular details. The blue sphere at the center of each glucose monomer is a formatting preference of mine, providing an easy way to identify the monosaccharide, a format known as Symbol Nomenclature for Glycans (Varki et al. 2015). If you would like to utilize this representation in your research, here is a link to install this plugin in VMD (https://glycam.org/docs/othertoolsservice/2016/06/03/3d-symbol-nomenclature-for-glycans-3d-sng/index.html). 

### Collapsing the GMML2 Model in Vacuum

In this project, I want to emphasize that—to my knowledge—there are no experimental data that directly resolve the 3D molecular structure of CTA plastics. As a result, several assumptions must be made in order to construct a physically reasonable starting model.

First, I assume that glucose residues adopt the <sup>4</sup>C<sub>1</sub> chair conformation, based on (a) experimental crystal structures of related carbohydrates in the Protein Data Bank and (b) quantum mechanical calculations identifying <sup>4</sup>C<sub>1</sub> as the lowest-energy pucker. For readers unfamiliar with pyranose ring puckering notation: the number or letter shown as a superscript is above the plane of the ring, while subscripts are below. The unmarked “C” denotes a chair conformation. So <sup>4</sup>C<sub>1</sub> indicates that carbon 4 is above and carbon 1 is below the ring plane.

Second, I assume that the polymer chains are randomly coiled and entangled, in line with the general consensus in polymer science regarding the morphology of amorphous thermoplastics.

The next step is to pack multiple copies of the GMML2-generated polymer into a simulation box. This creates an “unentangled” model of the bulk plastic. I use PACKMOL for this packing step. While we could, in theory, pack the as-generated linear chains, doing so would require a prohibitively large simulation box to prevent overlap between the extended polymers—leaving vast regions of vacuum. This would not only lead to long equilibration times, but also make chain entanglement unlikely or unphysical within practical simulation timescales.

Instead, a more effective approach is to first collapse each polymer chain in vacuum at high temperature, and then pack these compacted structures together to form the initial amorphous system.

Now, you might be thinking, “*It looks like all the Glc monomers are already in the <sup>4</sup>C<sub>1</sub> conformation*,” which is true—and seems very convenient at first glance. Unfortunately, it’s been our experience that bond angle force constants within pyranose rings in the GLYCAM06 force field are artificially weak, often resulting in small populations of unrealistic ring conformations even at ambient temperature. These artifacts are likely exacerbated at high temperature, where the polymers have greater thermal energy to overcome the barrier into nonphysical puckers.

Fortunately, this will be addressed in our upcoming force field release—GLYCAM25, coming soon!

For the time being, a workaround, albeit a bit tedious, is to apply NMR-style torsion angle restraints to all of the rings in the polymer to maintain the <sup>4</sup>C<sub>1</sub> pucker during collapse.

As an example, let's say you have a butane molecule and you want to restrain the C1-C2-C3-C4 torsion of this molecule to 180&deg; during an AMBER MD simulation. The restraint file would look like this:
```
&rst
iat=1,2,3,4,
r1=-160.0, r2=-170.0, r3=170.0, r4=160.0,
rk2=30.0, rk3=30.0,
&end
```
In your torsion restraint file, each restraint block must start with "$rst" and end with "&end". Unfortunately, restraint blocks can only be specified by atom number, which is why it's so tedious, and that is what the flag "iat" represents. I pray that someone will edit the source code in the future to allow for atom name masks instead of just atom number, but we will have to make do. Sander will apply the torsion restraint via a harmonic well potential, where r1 represents the most negative and r4 represents the most positive torsion defining your well. In practice, a solid well shape is such that each r is seperated by 10 degrees, with your angle bisecting r2 and r3. This is exactly what I have above to restrain the torsion to 180&deg;. 



For now, here is a short video of the collapse, with CTA shown in twister format in VMD:

<video width="640" height="360" controls>
  <source src="/figures/CTA_collapse.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

**In my next blog update, I will provide specifics on how to run the collapse in AMBER with NMR torsion restraints.**

*Last Updated 6/27/25*

## References
[1] Ajit Varki, Richard D. Cummings, Markus Aebi, Nicole H. Packer, Peter H. Seeberger, Jeffrey D. Esko, Pamela Stanley, Gerald Hart, Alan Darvill, Taroh Kinoshita, James J. Prestegard, Ronald L. Schnaar, Hudson H. Freeze, Jamey D. Marth, Carolyn R. Bertozzi, Marilynn E. Etzler, Martin Frank, Johannes F. G. Vliegenthart, Thomas Lütteke, Serge Perez, Evan Bolton, Pauline Rudd, James Paulson, Minoru Kanehisa, Philip Toukach, Kiyoko F. Aoki-Kinoshita, Anne Dell, Hisashi Narimatsu, William York, Naoyuki Taniguchi, and Stuart Kornfeld. 2015. Symbol nomenclature for graphical representations of glycans. Glycobiology 25, 12, https://doi.org/10.1093/glycob/cwv091.



