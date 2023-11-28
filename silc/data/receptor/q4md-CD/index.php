<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title> R.E.DD.B.: code F-85 </title>
<meta http-equiv="Content-type" content="text/html;charset=iso-8859-1" />
<link rel="stylesheet" media="screen" type="text/css" title="Essai" href="http://upjv.q4md-forcefieldtools.org/REDDB/up.css" />
<link rel="shortcut icon" href="http://upjv.q4md-forcefieldtools.org/REDDB/images/favicon.ico" type="image/x-icon" />
<script type="text/javascript" src="http://upjv.q4md-forcefieldtools.org/Jmol/Jmol.js"></script>
<script type="text/javascript" src="../../../jquery/jquery-1.3.2.min.js"></script>

<script type="text/javascript">
function GetId(id)
{
return document.getElementById(id);
}
var i=false; // La variable i nous dit si la bulle est visible ou non
      
function move(e)
{
if(i)
{  // Si la bulle est visible, on calcul en temps reel sa position ideale
if (navigator.appName!="Microsoft Internet Explorer")
{ // Si on est pas sous IE
GetId("curseur").style.left=e.pageX + 5+"px";
GetId("curseur").style.top=e.pageY + 10+"px";
}
else
{ // Modif propose par TeDeum, merci a lui
if(document.documentElement.clientWidth>0)
{
GetId("curseur").style.left=20+event.x+document.documentElement.scrollLeft+"px";
GetId("curseur").style.top=10+event.y+document.documentElement.scrollTop+"px";
} 
else
{
GetId("curseur").style.left=20+event.x+document.body.scrollLeft+"px";
GetId("curseur").style.top=10+event.y+document.body.scrollTop+"px";
}
}
}
}
function montre(text)
{
if(i==false) {
GetId("curseur").style.visibility="visible"; // Si il est cache (la verif n'est qu'une securite) on le rend visible.
GetId("curseur").innerHTML = text; // on copie notre texte dans l'element html
i=true;
}
}

function cache()
{
if(i==true)
{
GetId("curseur").style.visibility="hidden"; // Si la bulle est visible on la cache
i=false;
}
}

document.onmousemove=move; // des que la souris bouge, on appelle la fonction move pour mettre a jour la position de la bulle.
</script>

<script type="text/javascript">
function popup (page,nom,largeur,hauteur,options)
{
var top = (screen.height-hauteur)/2;
var left = (screen.width-largeur)/2;
window.open (page,nom,"top="+top+",left = "+left+",width="+largeur+",height="+hauteur+","+options);
}
</script>

</head>
<body>
<div id="curseur" class="infobulle"></div> <!-- bulle d information -->
<p>
<img src="../../images/REDDB2.gif" style="float: left;" alt="" /><img src="../../images/REDDB2.gif" style="float: right;" alt="" />
</p>

<div class="centre">
<span class="r5">
<br /><br />Summary of information<br />
PROJECT<br />Cyclodextrins<br /><br />
</span>
<p>
<span class="centre">
<a href="http://en.wikipedia.org/wiki/Cyclodextrins" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/Wikipedia.jpg" alt="Wiki Logo" height="100" /></a><br /><br />
</span>
</p>
</div>

<div class="contenu">
<p>
<span class="value3"><b>Upload</b></span><b> 24-11-2009 13:37 (Day-Month-Year, Paris time)</b><br />
<span class="value3"><b>Update</b></span><b> 27-09-2011 13:29 (Day-Month-Year, Paris time)</b><br />
</p>
<span class="bl2">
<br /><br /><b>Information about the Author (who submitted the project in R.E.DD.B.)</b><br /><br />
</span>
<p><span class="value"><b>Firstname</b></span> Christine</p>
<p><span class="value"><b>Lastname</b></span> Cezard</p>
<p><span class="value"><b>Institute</b></span> UFR de pharmacie, UPJV</p>
<p><span class="value"><b>City</b></span> Amiens</p>
<p><span class="value"><b>Country</b></span> FRANCE<br /><br /><br /></p>

<span class="bl2">
<b>General information about the Project</b><br /><br />
</span>
<div class="b2">
<p><b>Molecule keywords</b></p>

<table class="keywords">
<tr>
<td class="pair">Organoglycopeptide, Glycoconjugate</td>
<td class="impair">Cycloamylose</td>
<td class="pair">Oligosaccharide, Glucoside, Mannoside</td>
<td class="impair">Organic protective groups: OMe, OAc, OBn, OBz</td>
<td class="pair">Ligand carrier</td>
</tr>
</table>
</div>
<p><br /><b>Abstract</b></p>
<div class="summary">
<i>Introduction</i> <br />
Cyclodextrins (CDs) are macrocyclic oligosaccharides composed of <i>D</i>-glucose units linked together by &#945;(1&#8594;4) glycosidic bonds. CDs are versatile compounds with various applications in chemistry and biology: they can be used as building blocks for supramolecular structures and/or act as drug/ligand carriers as one of their ability is to stabilize soluble organic molecules in liquid phases. To affect their solubility and/or their selectivity towards specific targets, CDs are substituted. Derivatization of native CDs can be of several types: (i) substitution with common organic protecting groups to modify the solubility properties and (ii) grafting of ligands such as amino-acid derived arms to modulate recognition mechanisms. In the absence within the GLYCAM force field of the fragments required to build the aforementioned substituted CDs, a new Force Field Topology DataBase (FFTopDB) compatible with the reported CD based molecular systems has been developed. In this work, the considered substitutents are heterogeneous, <i>i. e.</i> of organic, peptidic and carbohydrate nature.<br />
<div style="text-align:center;"><img src="./script4.ff" alt="" /><br /><span style="font-size:12px;">Scheme 1</span><br /></div>
<span style="font-size:12px;">Multiple orientation, multiple conformation and multiple molecule charge derivation and force field library building for cyclodextrin based systems automatically carried out using the R.E.D. IV program. (A) Description of the eight molecules involved in the procedure and construction of the corresponding FFTopDB. (Dashed) light-grey lines: intra-molecular charge constraints used during the fit; (dashed) dark-grey lines: inter-molecular charge constraints used during the fit. (B) Building of native and substituted CDs using the FFtopDB generated.<br /><br /></span>
<i>Computational details</i><br />
Multiple molecules, multiple conformations and multiple molecular orientations were used during the charge calculation procedure. RESP charge derivation and force field library building for the presented fragments were carried out in a single R.E.D. IV job following the procedure summarized on Scheme 1 and corresponding to a 34-structure RESP fit. Eight molecules, namely methyl &#945;-<i>D</i>-glucopyranoside, &#946;-<i>N</i>-acetamido-<i>D</i>-mannoside, methyl &#945;-<i>D</i>-mannoside, dimethyl ether, methyl acetate, benzyl methyl ether, methyl benzoate and <i>N,N&#39;</i>-dimethylsuccinamide were considered. Two conformations were selected and optimized for the sugar-derived units since Glucose and Mannose monosaccharide derivatives show two main populations in solution relative to a <i>gg/gt</i> isomerism around the the &#969; dihedral angle. Four molecular orientations based on the (C1 C3 C5) and (C2 C4 O5) sets of atoms were considered in molecular electrostatic potential (MEP) computation.<a href="#ref1">[1]</a> For each of the protecting groups/ligands, only one conformation corresponding to the lowest minimum obtained after geometry optimization was taken into account. Two orientations, based on the connected C(Me)-C-O atoms for the organic protecting groups and based on the (C1 C2 C3) atoms for <i>N,N&#39;</i>-dimethylsuccinamide were considered.<a href="#ref1">[1]</a> No weighting factor to bias a conformation over another was applied during the charge derivation procedure. Inter- and intra-molecular charge constraints were used during the fitting step to define the required molecular fragments. Inter-molecular charge constraints were set to a target value of zero between (i) the 2-hydroxyl, 3-hydroxyl and 6-hydroxyl groups of the methyl &#945;-<i>D</i>-glucopyranoside and the methoxy group belonging to the different organic protecting groups and (ii) the 3-hydroxyl, 4-hydroxyl and 6-hydroxyl groups of the &#946;-<i>N</i>-acetamido-<i>D</i>-mannoside and the methyl group of methyl &#945;-<i>D</i>-mannoside. Intra-molecular charge constraints between the 4-hydroxyl and the methyl group of the methyl &#945;-<i>D</i>-glucopyranoside as well as for one of the <i>NH</i>-methyl group of <i>N,N&#39;</i>-dimethylsuccinamide and the acetyl group connected to &#946;-<i>N</i>-acetamido-<i>D</i>-mannoside were set to zero. For the sake of compatibility between organic/sugar and peptidic/sugar connections, an additional intra-molecular charge constraint set to a value of 0.1980 was imposed for the methyl group of <i>N,N&#39;</i>-dimethylsuccinamide (that not already involved in the intra-molecular charge constraint previously defined). Geometry optimization and MEP computation were carried out with the Gaussian03 program, while charge fitting was done using the RESP program. Geometry optimization was carried out at the RHF/6-31G** level of theory, whereas MEP computation was performed at the RHF/6-31G* one using the Connolly surface algorithm. The molecular orientation of optimized geometries was controlled before MEP calculation using the rigid-body reorientation algorithm implemented in R.E.D. leading to highly reproducible charge values. RESP charge fitting was carried out following a two-stage fitting procedure with a hyperbolic restraint function, using a weighting factor of 0.0005 and 0.001 for the two stages, respectively. A RRMS (Relative Root Mean Square) value of 0.096 was obtained for the charge fitting step between the MEP calculated by quantum chemistry and that generated using the derived charge values. <br /><br />
<i>Discussion</i> <br />
CDs substituted with organic protecting groups or with an amino-acid derived arm can be considered like heterogeneous systems as they display carbohydrate, organic and peptidic components. The modeling of such a molecular system should be heterogeneous as well and different force fields specific to each of these constituents are commonly used. The GLYCAM approach for deriving atomic charges has proven to be highly innovatrice and effective: the weighting of the different conformations taken into account during the charge fitting step is related to their occurence during molecular dynamics (MD) simulations.<a href="#ref2">[2]</a> However, this approach applied to study the conformational space of a ligand to be grafted on a CD remains a complex task. Consequently, considering that the highly heterogeneous CDs studied in this work are of a &#147;glyco-organo-peptidic&#148; nature, we have chosen to derive RESP charges following the &#147;Amber&#148; approach as a single and homogeneous procedure. The charge derivation and force field library building procedures only involve representative minima optimized by quantum mechanics in the gas phase followed by a MEP computation step using the Connolly surface algorithm and a two RESP-charge fitting stage.<a href="#ref3">[3]</a> Furthermore, for the sake of consistency and homogeneity, we developed a unique force field, namely &#147;q4md-CD&#148;, to model these CD based systems. q4md-CD presents the following features: (i) atom charges for all required fragments were derived following the procedure previously reported. The glycocluster FFTopDB is available as a suite of force field libraries for molecular fragments in the Tripos mol2 file format and is used to build Amber OFF force field libraries. A LEaP <a href="./script1.ff" onclick="window.open(this.href,'_blank');return false;">script</a> is available for this purpose in this project. (ii) A minimal number of charge constraints was used during the fitting step to keep the bias on the RRMS as low as possible, and consequently charge values of aliphatic hydrogen atoms were not constrained to the zero value. (iii) Scaling factor values of 1.2 and 2.0 for the 1-4 electrostatic and van der Waals interactions were used in MD simulations, respectively. (iv) Geometrical parameters are adapted to the nature of the considered ligand. Force field parameters used in association with the FFTopDB reported are described in the <a href="./script3.ff" onclick="window.open(this.href,'_blank');return false;">frcmod</a> file. (v) Compatibility with amino acid residues from the Amber FFTopDB is consequently assured.<br /><br />
<i>Conclusion</i> <br />
Derived charge values have been validated after analyses of 50 nsec MD simulations. Structural characteristics (hydrogen bond patterns, &#969; dihedral populations, sugar pucker and CD solvation) observed during MD compare well to experimental data. To the best of our knowledge q4md-cd represents the first development of an Amber-GLYCAM hybrid force field, which allows studying glycoconjugates in a fully homogeneous approach. Extension of the force field to more complex glycoconjugates is underway. The corresponding study will be included in a forthcoming paper.<br /><br />
<span style="font-size:12px;"><a name="ref1">[1]</a>&nbsp; Molecular orientations for monosaccharide units are based on the (C1 C3 C5), (C5 C3 C1), (C2 C4 O5) and (O5 C4 C2) sets of three atoms. Molecular orientations for the organic protecting groups are based on the: (i) (CX OS CM) and (CM OS CX) sets of three atoms for dimethyl ether, (ii) (CX OS C) and (C OS CX) sets of three atoms for methyl acetate, (iii) (CX OS CM) and (CM OS CX) sets of three atoms for  benzyl methyl ether, and (iv) (CX OS C) and (C OS CX) sets of three atoms for methyl benzoate. Molecular orientations for <i>N,N&#39;</i>-dimethylsuccinamide are based on the (C1 C2 C3) and (C3 C2 C1) sets of three atoms. See the corresponding PDB files for atom naming convention.<br />
<a name="ref2">[2]</a>&nbsp; Woods, R. J.; Chappelle, R. <a href="http://www.sciencedirect.com/science?_ob=ArticleURL&amp;_udi=B6TGT-40WDT97-J&amp;_user=781134&amp;_coverDate=08%2F04%2F2000&amp;_rdoc=16&amp;_fmt=high&amp;_orig=browse&amp;_srch=doc-info(%23toc%235263%232000%23994729998%23205401%23FLA%23display%23Volume)&amp;_cdi=5263&amp;_sort=d&amp;_docanchor=&amp;_ct=28&amp;_acct=C000043238&amp;_version=1&amp;_urlVersion=0&amp;_userid=781134&amp;md5=d8df08b1b2188474ecf5d7df18c9f10f" onclick="window.open(this.href,'_blank');return false;"><i>J. Mol. Struct. (Theochem)</i> 2000, 527, 149-156</a>.&nbsp; Basma, M.; Sundara, S.; Calgan, D.; Varnali, T.; Woods R. J. <a href="http://www3.interscience.wiley.com/journal/82004191/abstract" onclick="window.open(this.href,'_blank');return false;"><i>J. Comput. Chem.</i> 2001, 22, 1125-1137</a>.&nbsp; Kirschner, K. N.; Yongye, A. B.; Tschampel, S. M.; Gonzalez-Outeirino, J.; Daniels, C. R.; Lachele Foley, B.; Woods, R. J. <a href="http://www3.interscience.wiley.com/journal/116311717/abstract" onclick="window.open(this.href,'_blank');return false;"><i>J. Comput. Chem.</i> 2007, 29, 622-655</a>.<br />
<a name="ref3">[3]</a>&nbsp; Bayly, C. I.; Cieplak, P.; Cornell, W. D.; Kollman, P. A. <a href="http://pubs.acs.org/doi/abs/10.1021/j100142a004" onclick="window.open(this.href,'_blank');return false;"><i>J. Phys. Chem.</i> 1993, 97, 10269-10280</a>.&nbsp; Cornell, W. D.; Cieplak, P.; Bayly, C. I.; Kollman, P. A. <a href="http://pubs.acs.org/doi/abs/10.1021/ja00074a030?journalCode=jacsat&amp;quickLinkVolume=115&amp;quickLinkPage=9620&amp;volume=115" onclick="window.open(this.href,'_blank');return false;"><i>J. Am. Chem. Soc.</i> 1993, 115, 9620-9631</a>.&nbsp; Cieplak, P.; Cornell, W. D.; Bayly, C. I.; Kollman, P. A. <a href="http://www3.interscience.wiley.com/journal/109583237/abstract" onclick="window.open(this.href,'_blank');return false;"><i>J. Comput. Chem.</i> 1995, 16, 1357-1377</a>.</span>

</div>
<br />
<p><span class="value"><b>Publication</b></span> YES &nbsp; &nbsp; &nbsp; <a href="http://www.ncbi.nlm.nih.gov/pubmed/21792425?dopt=Abstract" onclick="window.open(this.href,'_blank');return false;"><img src="../../images/pubmed.png" alt="" /></a></p>
<p><span class="value"><b>Author(s)</b></span> C. Cezard, X. Trivelli, F. Aubry, F. Djedaini-Pilard and F.-Y. Dupradeau</p>
<p><span class="value"><b>Journal</b></span><i> Phys. Chem. Chem. Phys.</i></p>
<p><span class="value"><b>Year</b></span><b> 2011</b></p>
<p><span class="value"><b>Volume</b></span><i> 13</i></p>
<p><span class="value"><b>Page(s)</b></span> 15103-15121</p><br />
<span class="value2"><b>"Whole molecule" or "Molecule fragment" type project</b></span> MOLECULE FRAGMENT <br /><br />
<span class="value2"><b>Interface R.E.D. used ?</b></span> YES<br />
<span class="bl2">
<br /><br /><b>Charge derivation procedure</b><br /><br />
</span>
<p>
<span class="value2"><b>Number of Tripos mol2 file(s) provided by the author(s)</b></span><span class="flash"> 22</span><br />
</p>
<p>
Contain charge values &#38; information about molecular topology<br /><br />
</p>

<div class="information">
<table id="info">
<thead>
<tr>
<th onmouseover="montre('Molecule or fragment number');" onmouseout="cache();">No</th>
<th onmouseover="montre('Name of each molecule or fragment');" onmouseout="cache();">Name</th>
<th onmouseover="montre('Download each force field library individually');" onmouseout="cache();">Download</th>
<th onmouseover="montre('Get general information about each molecule or fragment');" onmouseout="cache();">Wikipedia</th>
<th onmouseover="montre('Use Jmol and JRE(sun) to display each structure');" onmouseout="cache();">3D visualization</th>
</tr>
</thead>
<tbody>
<tr class="impair">
<td>1</td>
<td>Fragment MG2</td>
<td><a href="tripos1.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MG2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-1.php","JavaApplet1",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>2</td>
<td>Fragment MG3</td>
<td><a href="tripos2.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MG3" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-2.php","JavaApplet2",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>3</td>
<td>Fragment MG6</td>
<td><a href="tripos3.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MG6" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-3.php","JavaApplet3",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>4</td>
<td>Fragment MGA</td>
<td><a href="tripos4.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MGA" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-4.php","JavaApplet4",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>5</td>
<td>Fragment MGB</td>
<td><a href="tripos5.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MGB" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-5.php","JavaApplet5",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>6</td>
<td>Fragment MGC</td>
<td><a href="tripos6.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MGC" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-6.php","JavaApplet6",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>7</td>
<td>Fragment MGO</td>
<td><a href="tripos7.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MGO" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-7.php","JavaApplet7",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>8</td>
<td>Fragment MGR</td>
<td><a href="tripos8.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment MGR" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-8.php","JavaApplet8",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>9</td>
<td>Fragment OAC</td>
<td><a href="tripos9.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment OAC" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-9.php","JavaApplet9",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>10</td>
<td>Fragment OBN</td>
<td><a href="tripos10.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment OBN" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-10.php","JavaApplet10",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>11</td>
<td>Fragment OBZ</td>
<td><a href="tripos11.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment OBZ" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-11.php","JavaApplet11",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>12</td>
<td>Fragment OME</td>
<td><a href="tripos12.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment OME" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-12.php","JavaApplet12",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>13</td>
<td>Fragment SCC</td>
<td><a href="tripos13.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment SCC" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-13.php","JavaApplet13",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>14</td>
<td>Fragment BM3</td>
<td><a href="tripos14.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BM3" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-14.php","JavaApplet14",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>15</td>
<td>Fragment BM4</td>
<td><a href="tripos15.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BM4" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-15.php","JavaApplet15",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>16</td>
<td>Fragment BM6</td>
<td><a href="tripos16.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BM6" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-16.php","JavaApplet16",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>17</td>
<td>Fragment BMA</td>
<td><a href="tripos17.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BMA" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-17.php","JavaApplet17",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>18</td>
<td>Fragment BMB</td>
<td><a href="tripos18.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BMB" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-18.php","JavaApplet18",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>19</td>
<td>Fragment BMC</td>
<td><a href="tripos19.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BMC" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-19.php","JavaApplet19",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>20</td>
<td>Fragment BMO</td>
<td><a href="tripos20.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BMO" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-20.php","JavaApplet20",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="impair">
<td>21</td>
<td>Fragment BMR</td>
<td><a href="tripos21.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment BMR" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-21.php","JavaApplet21",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
<tr class="pair">
<td>22</td>
<td>Fragment AMO</td>
<td><a href="tripos22.mol2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Fragment AMO" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
<td><a href='javascript:popup("http://upjv.q4md-forcefieldtools.org/REDDB/projects/F-85/JavaApplet-22.php","JavaApplet22",750,500,"menubar=no,scrollbars=no,statusbar=no,resizable=yes")'><img style="margin-left: auto;margin-right: auto;" src="../../images/Jmol.png" alt="Jmol Logo" /></a></td>
</tr>
</tbody>
</table><br /><br />  
<span class="value2"><b>Number of molecule(s) used in the charge derivation procedure</b></span><span class="flash"> 8</span>
<p>
File(s) provided to the PDB format<br /><br />
</p>
<table class="info">
<thead>
<tr>
<th onmouseover="montre('Molecule number');" onmouseout="cache();">No</th>
<th onmouseover="montre('Name of each molecule used in charge derivation');" onmouseout="cache();">Molecule name</th>
<th onmouseover="montre('Number of conformation(s) for each molecule used in charge derivation');" onmouseout="cache();">Conformation No</th>
<th onmouseover="montre('Re-orientation procedure applied before MEP computation for each conformation');" onmouseout="cache();">Reorientation procedure</th>
<th onmouseover="montre('Number of molecular orientation(s) used in MEP computation for each conformation');" onmouseout="cache();">Mol. orientation No</th>
<th onmouseover="montre('Download each molecule used in charge derivation individually');" onmouseout="cache();">Download</th>
<th onmouseover="montre('Get general information about each molecule');" onmouseout="cache();">Wikipedia</th>
</tr>
</thead>
<tbody>
<tr class="impair">
<td>1</td>
<td>Methyl alpha-D-glucopyranoside</td>
<td>2</td>
<td>Rigid Body Reorient Algo</td>
<td>4</td>
<td><a href="mol1.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Methyl alpha-D-glucopyranoside" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="pair">
<td>2</td>
<td>Methyl acetate</td>
<td>1</td>
<td>Rigid Body Reorient Algo</td>
<td>2</td>
<td><a href="mol2.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Methyl acetate" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="impair">
<td>3</td>
<td>Benzyl methyl ether</td>
<td>1</td>
<td>Rigid Body Reorient Algo</td>
<td>2</td>
<td><a href="mol3.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Benzyl methyl ether" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="pair">
<td>4</td>
<td>Methyl benzoate</td>
<td>1</td>
<td>Rigid Body Reorient Algo</td>
<td>2</td>
<td><a href="mol4.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Methyl benzoate" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="impair">
<td>5</td>
<td>Dimethyl ether</td>
<td>1</td>
<td>Rigid Body Reorient Algo</td>
<td>2</td>
<td><a href="mol5.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Dimethyl ether" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="pair">
<td>6</td>
<td>N,N-dimethylsuccinamide</td>
<td>1</td>
<td>Rigid Body Reorient Algo</td>
<td>2</td>
<td><a href="mol6.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/N,N-dimethylsuccinamide" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="impair">
<td>7</td>
<td>beta-N-acetamido-D-mannoside</td>
<td>2</td>
<td>Rigid Body Reorient Algo</td>
<td>4</td>
<td><a href="mol7.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/beta-N-acetamido-D-mannoside" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
<tr class="pair">
<td>8</td>
<td>Methyl alpha-D-mannoside</td>
<td>2</td>
<td>Rigid Body Reorient Algo</td>
<td>4</td>
<td><a href="mol8.pdb" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></td>
<td><a href="http://en.wikipedia.org/wiki/Methyl alpha-D-mannoside" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/wiki.png" alt="Wiki Logo" /></a></td>
</tr>
</tbody>
</table>
</div>

<p>
<span class="bl2">
<br /><br /><b>Information regarding Quantum Calculations</b><br /><br />
</span>
</p>
<p>
<b><i>Geometry optimization</i></b>
</p>
<p><span class="value"><b>Program 1</b></span> GAUSSIAN 2003</p>
<p><span class="value"><b>Theory level 1</b></span> HF </p>
<p><span class="value"><b>More information 1</b></span> Opt=Tight</p>
<p><span class="value"><b>Basis set 1</b></span> 6-31G**<br /><br /></p>
<p>
<b><i>Molecular electrostatic potential computation</i></b>
</p>
<p><span class="value"><b>Program 2</b></span> GAUSSIAN 2003</p>
<p><span class="value"><b>Theory level 2</b></span> HF </p>
<p><span class="value"><b>More information 2</b></span> IOp(6/33=2) NoSymm</p>
<p><span class="value"><b>Basis set 2</b></span> 6-31G*</p>
<p><span class="value"><b>Algorithm</b></span> CONNOLLY SURFACE<br /><br /><br /></p>
<span class="bl2">
<b>Information about the charge fit</b><br /><br />
</span>
<p><span class="value"><b>Program</b></span> RESP</p>
<p><span class="value"><b>Number of stage(s)</b></span> 2</p>
<p><span class="value">input of stage 1 </span><a href="input1.in" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></p>
<p><span class="value">input of stage 2 </span><a href="input2.in" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a></p><br /><span class="bl2">
<br /><b>Files the author of the project wishes to provide...</b><br /><br />
</span>
<span class="value2">A script to convert Tripos mol2 file(s) into LEaP OFF library(ies) (for AMBER)...</span><a href="./script1.ff" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br />
<span class="value2">A script to convert Tripos mol2 file(s) into RTF or PSF library(ies) (for CHARMM)...</span><a href="./script2.ff" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br />
<span class="value2">A file to provide new force field parameters compatible with the Tripos mol2 file(s)...</span><a href="./script3.ff" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br />
<span class="value2">A file (choice made by the author) to provide more information about the project...</span><a href="./script4.ff" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br />
<span class="value2">A file (choice made by the author) to provide more information about the project...</span><a href="./script5.ff" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br /><br />

<span class="bl2">
<span class="value2"><b>Download the whole project...</b></span>
</span>
<a href="./F-85.tar.bz2" onclick="window.open(this.href,'_blank');return false;"><img style="margin-left: auto;margin-right: auto;" src="../../images/disquette.gif" alt="Link" /></a><br /><br /><br /><br />

<p style="text-align: center;">
<a href="http://validator.w3.org/check?uri=referer" onclick="window.open(this.href,'_blank');return false;"><img src="http://www.w3.org/Icons/valid-xhtml10" alt="Valid XHTML 1.0 Strict" height="31" width="88" /></a>
<a href="http://jigsaw.w3.org/css-validator/check/referer" onclick="window.open(this.href,'_blank');return false;"><img style="border:0;width:88px;height:31px" src="http://jigsaw.w3.org/css-validator/images/vcss" alt="CSS Valide !" /></a><br /><br />
<span class="spec2">
<?php include('../../official-date.php'); ?> <?php include('../../official-institute.php'); ?>R.E.DD.B. projects <a href='http://www.gnu.org/philosophy/philosophy.html' onclick="window.open(this.href,'_blank');return false;">free</a>.<br /><br />
</span>
</p>
</div>
</body>
</html>
