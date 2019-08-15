# PROJ-DATUMGRID-EUROPE

The files in this package can be unpacked in the PROJ data directory.
For an installed PROJ this may be /usr/local/share/proj or /usr/share/proj
on unix style operating systems.

## About the Europe PROJ resource package

The Europe PROJ resource package is a collections of grids
and init-files that are not essential for the functionality of PROJ but
still of general interest to PROJ-users in the region specific to the
package.

## Included grids

### Denmark: DVR90 heights -> ETRS89 ellipsoidal heights

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transform local height system DVR90 (EPSG:5799) to ETRS89 ellipsoidal heights (EPSG:4937). The horizontal
grid coordinates are referenced to ETRS89.

* dvr90.gtx

### Denmark: DNN heights -> ETRS89 ellipsoidal heights

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transform local height system DNN (EPSG:5733) to ETRS89 ellipsoidal heights (EPSG:4937). The horizontal
grid coordinates are referenced to ETRS89.

* dnn.gtx

### Faroe Islands: FVR09 heights -> ETRS89 ellipsoidal heights

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transform local height system FVR09 (EPSG:5317) to ETRS89 ellipsoidal heights (EPSG:4937). The horizontal
grid coordinates are referenced to ETRS89.

* fvr09.gtx

### Nordic+Baltic countries: Deformation model

*Source*: [The Nordic Geodetic Commission](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: CTable2 and GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Deformation model covering the Nordic and Baltic countries. Used in
transformations between global reference frames and the local realisations
of ETRS89 in the Nordic and Baltic countries. See the `NKG` init-file for
examples of use.

* nkgrf03vel_realigned_xy.ct2
* nkgrf03vel_realigned_z.gtx

### Ireland: OSGM15 height, Malin head datum -> ETRS89 ellipsoidal heights

*Source*: [Ordnance Survey](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html)  
*Format*: GTX  
*License*: [The 2-Clause BSD License](https://opensource.org/licenses/bsd-license.php)  

Vertical transformation for Geoid model OSGM15, Malin head datum (EPSG:5731). Used in
transformation from OSGM15 orthometric heights to WGS84/ETRS89 ellipsoidal heights.
The Malin Head datum is used in the Republic of Ireland.

* OSGM15_Malin.gtx

### Northern Ireland: OSGM15 height, Belfast height -> ETRS89 ellipsoidal heights

*Source*: [Ordnance Survey](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html)  
*Format*: GTX  
*License*: [The 2-Clause BSD License](https://opensource.org/licenses/bsd-license.php)  

Vertical transformation for Geoid model OSGM15, Belfast datum (EPSG:5732). Used in
transformation from OSGM15 orthometric heights to WGS84/ETRS89 ellipsoidal heights.
The Belfast datum is used in the Northern Ireland.

* OSGM15_Belfast.gtx

### United Kingdom: OSGB36 -> ETRS89

*Source*: [Ordnance Survey](https://www.ordnancesurvey.co.uk/business-and-government/help-and-support/navigation-technology/os-net/formats-for-developers.html)  
*Format*: NTv2  
*License*: [The 2-Clause BSD License](https://opensource.org/licenses/bsd-license.php)

OSTN15 is the definitive OSGB36/ETRS89 transformation. OSTN15 in combination with the ETRS89
coordinates of the OS Net stations, rather than the old triangulation network, define the National
Grid. This means that, for example, the National Grid coordinates of an existing OSGB36 point,
refixed using GNSS from OS Net and OSTN15, will be the correct ones. The original archived
OSGB36 National Grid coordinates of the point (if different) will be wrong, by definition, but the
two coordinates (new and archived) will agree on average to better than 0.1m (0.1m rmse, 68%
probability).

* OSTN15_NTv2_OSGBtoETRS.gsb

### Belgium: BD72 -> ETRS89

Grid transformation from Belgium Datum 72 to ETRS89 in Belgium.

*Source*: [IGN](http://www.ngi.be/Common/Lambert2008/NTv2.zip)
*Format*: NTv2
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)
*Credit*: (c)2014-2016 Grid created by Nicolas SIMON.

It provides an accuracy better than 5mm in 99.6% of validation points and a
worst case of 1.3 cm outside the border.
Recommended to transform from Belgian Lambert 72 (EPSG:31370) to Belgian Lambert 2008 (EPSG:3812)
Documentation in French: (http://www.ngi.be/FR/FR2-1-7.shtm)
Documentation in Dutch:  (http://www.ngi.be/NL/NL2-1-7.shtm)

* bd72lb72_etrs89lb08.gsb

### French vertical grids

*Source*: [IGN France](https://geodesie.ign.fr/index.php?page=grilles)  
*Converter*: build_french_vgrids.sh  
*Format*: GTX  
*License*: [License Ouverte (in French)](https://www.etalab.gouv.fr/wp-content/uploads/2017/04/ETALAB-Licence-Ouverte-v2.0.pdf) / [Open License (English translation)](https://www.etalab.gouv.fr/wp-content/uploads/2014/05/Open_Licence.pdf)

Grids to convert from national vertical datum to ellipsoidal height.

Continental France an Corsica:
* RAF09.gtx: Continental France, NGF-IGN 1969 to RGF93 (ETRS89) (2009 version)
* RAF18.gtx: Continental France, NGF-IGN 1969 to RGF93 (ETRS89) (2018 version)
* RAC09.gtx: Corsica, NGF-IGN 1978 to RGF93 (ETRS89)

French Antilles to RGAF09 (Réseau Géodésique des Antilles Françaises 2009):
* RAGTBT2016.gtx: Guadeloupe, IGN 1988 (GUADELOUPE) to RGAF09
* RALD2016.gtx: La Desirade, IGN 2008 LD (GUADELOUPE / LA DESIRADE) to RGAF09
* RALS2016.gtx: Les Saintes, IGN 1988 LS (GUADELOUPE / LES SAINTES) to RGAF09
* RAMART2016.gtx: Martinique, IGN 1987 (MARTINIQUE) to RGAF09
* RAMG2016.gtx: Marie Galante, IGN 1988 MG (GUADELOUPE / MARIE-GALANTE) to RGAF09
* gg10_sbv2.gtx: Saint Barthelemy, IGN 1988 SB (GUADELOUPE / SAINT-BARTHELEMY) to RGAF09
* gg10_smv2.gtx: Saint Martin, IGN 1988 SM (GUADELOUPE / SAINT-MARTIN) to RGAF09

French Antilles to WGS 84 (RRAF: Réseau de Référence des Antilles Françaises):
* ggg00v2.gtx: Guadeloupe, IGN 1988 (GUADELOUPE) to WGS 84 (RRAF)
* RALDW842016.gtx: La Desirade, IGN 2008 LD (GUADELOUPE / LA DESIRADE) to WGS 84 (RRAF)
* ggg00_lsv2.gtx: Les Saintes, IGN 1988 LS (GUADELOUPE / LES SAINTES) to WGS 84 (RRAF)
* ggg00_mgv2.gtx: Marie Galante, IGN 1988 MG (GUADELOUPE / MARIE-GALANTE) to WGS 84 (RRAF)
* ggm00v2.gtx: Martinique, IGN 1987 (MARTINIQUE) to WGS 84 (RRAF)
* ggg00_sbv2.gtx: Saint Barthelemy, IGN 1988 SB (GUADELOUPE / SAINT-BARTHELEMY) to WGS 84 (RRAF)
* ggg00_smv2.gtx: Saint Martin, IGN 1988 SM (GUADELOUPE / SAINT-MARTIN) to WGS 84 (RRAF)

Other:
* ggguy15.gtx: Guyane, NIVELLEMENT GENERAL DE GUYANE (NGG) 1977 to RGFG95 (RESEAU GEODESIQUE FRANCAIS DE GUYANE 1995)
* ggker08v2.gtx: Iles Kerguelen, IGN 1962 (KERGUELEN) to RGTAAF07 (RESEAU GEODESIQUE DES TAAF 2007)
* ggm04v1.gtx: Mayotte GGM04, SHOM 1953 (MAYOTTE) to RGM04 (RESEAU GEODESIQUE DE MAYOTTE 2004)
* ggspm06v1.gtx: Saint-Pierre et Miquelon (GGSPM06), DANGER 1950 (SAINT-PIERRE-ET-MIQUELON) to RGSPM06 (RESEAU GEODESIQUE DE SAINT-PIERRE-ET-MIQUELON 2006)
* RAR07_bl.gtx: La Reunion, IGN 1989 (REUNION) to RGR92 (RESEAU GEODESIQUE DE LA REUNION 1992)

French Polynesia to RGPF (Réseau Géodésique de Polynésie Française)
* ggpf02-Bora.gtx: Bora bora, BORA_SAU 2001 to RGPF
* ggpf02-Huahine.gtx: Huahine, HUAHINE_SAU 2001 to RGPF
* ggpf02-Maiao.gtx: Maiao, MAIAO 2001 to RGPF
* ggpf02-Maupiti.gtx: Maupiti, MAUPITI_SAU 2001 to RGPF
* ggpf02-Raiatea.gtx: Raietea, RAIATEA_SAU 2001 to RGPF
* ggpf02-Tahaa.gtx: Tahaa, TAHAA_SAU 2001 to RGPF
* ggpf02-Tupai.gtx: Tupai, TUPAI 2001 to RGPF
* ggpf05-HivaOa.gtx: Hiva Oa, HIVA OA to RGPF
* ggpf05-Nuku.gtx: Nuku Hiva, NUKU HIVA ALTI to RGPF
* ggpf08-Fakarava.gtx: Fakarava, IGN 1966 (TAHITI) TO RGPF
* ggpf08-Gambier.gtx: Gambier, GAMBIER to RGPF
* ggpf08-Hao.gtx: Hao, HAO to RGPF
* ggpf08-Mataiva.gtx: Mataiva, MATAIVA to RGPF
* ggpf08-Raivavae.gtx: Raivavae, Raivavae to RGPF
* ggpf08-Reao.gtx: Reao, Reao to RGPF
* ggpf08-Rurutu.gtx: Rurutu, Ruruto to RGPF
* ggpf08-Tikehau.gtx: Tikehau, Tikehau to RGPF
* ggpf08-Tubuai.gtx: Tubuai, Tubuai to RGPF
* ggpf10-Moorea.gtx: Moorea, MOOREA 1981 (MOOREA_SAU 2001) to RGPF
* ggpf10-Tahiti.gtx: Tahiti, IGN 1966 (TAHITI) to RGPF

### Sweden: RH 2000 heights -> SWEREF 99 ellipsoidal heights

*Source*: [Lantmäteriet](https://www.lantmateriet.se/en/maps-and-geographic-information/GPS-och-geodetisk-matning/Om-geodesi/Transformationer/Geoidmodeller/)  
*Format*: GTX  
*License*: [Creative commons, CC0](https://creativecommons.org/publicdomain/zero/1.0/)

Vertical transformation for Geoid model SWEN17_RH2000. Used to make the
transitions from heights above sea level in RH 2000 (EPSG:5613) to heights above the
ellipsoid in SWEREF 99 (EPSG:4377). The RH 2000 height system is used on Swedish mainland, Öland
and Gotland islands.

* SWEN17_RH2000.gtx

## Included init-files

### Denmark: System definitions related to ETRS89

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Init files with various transformations and conversions of relevance in
Denmark. All definitions in the init file uses the local ETRS89 reference
frame as a pivot datum. More information can be found in the
[NordicTransformations](https://github.com/NordicGeodesy/NordicTransformations)
repository. While the `.pol` files are init-files they are not supposed to be
by themselves. The exist purely as helper files for the `DK` init-file.
The `.pol` files holds polynomial constants used for transformations between
Danish legacy coordinate reference systems and ETRS89.

* DK
* DK_bornholm.pol
* DK_bridges.pol
* DK_general.pol
* DK_jutland.pol
* DK_zealand.pol

### The Faroe Islands: System definitions related to ETRS89

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Init file with various transformations and conversions of relevance in
The Faroe Islands. All definitions in the init file uses the local ETRS89
reference frame as a pivot datum. More information can be found in the
[NordicTransformations](https://github.com/NordicGeodesy/NordicTransformations)
repository.

* FO

### Nordic+Baltic countries: ITRFxx -> Local ETRS89 realizations

*Source*: [The Nordic Geodetic Commission](https://github.com/NordicGeodesy/NordicTransformations)  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transformations to and from the common Nordic referenc frame NKG_ETRF00.
This init-file describes transformations between global reference frames
and NKG_ETRF00 as well as transformations from NKG_ERTF00 to the local
realisations of ETRS89 in each of the countries involved with NKG.

All transformations in this init-file uses the common Nordic frame as a
pivot datum. Exempt from this dogma are transformations with labels
starting with ``_``. Those transformations are "private" to this file and are
only used as steps in more complicated transformations.

* NKG
# PROJ-DATUMGRID-NORTH-AMERICA

The files in this package can be unpacked in the PROJ data directory.
For an installed PROJ this may be /usr/local/share/proj or /usr/share/proj
on unix style operating systems.

## About the North America PROJ resource package

The North America PROJ resource package is a collections of grids
and init-files that are not essential for the functionality of PROJ but
still of general interest to PROJ-users in the region specific to the
package.

## Included grids

### Canada: NAD27 -> NAD83

*Source*: [Natural Resources Canada](https://open.canada.ca/data/en/dataset/b3534942-31ea-59cf-bcc3-f8dc4875081a)  
*Format*: NTv2  
*License*: [Open Government Licence - Canada](http://open.canada.ca/en/open-government-licence-canada)

Transform between NAD27 and NAD83 in Canada.

* ntv2_0.gsb

### Greenland: GVR2016 height -> GR96 ellipsoidal heights

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transform local height system GVR2016 (EPSG:8267) to GR96 ellipsoidal heights (EPSG:4909). The horizontal
grid coordinates are referenced to GR96.

* gvr2016.gtx

### Greenland: GVR2000 height -> GR96 ellipsoidal heights

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*Format*: GTX  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Transform local height system GVR2000 (EPSG:8266) to GR96 ellipsoidal heights (EPSG:4909). The horizontal
grid coordinates are referenced to GR96.

* gvr2000.gtx

### USA: NAD83 -> NAD83 (HARN/HPGN)

*Source*: [.los/.las files coming from NOAA](https://www.ngs.noaa.gov/PC_PROD/NADCON/NADCON.zip)
*Format*: NTv2
*License*: Public Domain

Grid data for high precision conversion of geographic coordinates from
the original NAD83 (1986) (EPSG:4269) to NAD83 High Precision Geodetic Networks
(HPGNs) (EPSG:4152), also referred to as High Accuracy Reference Networks (HARNs).
NAD 83 coordinates based on the HPGN/HARN surveys changed approximately 0.2
to 1.0 meter relative to the original NAD 83 (1986) adjustment.

Those files have been converted to NTv2 grids with the
[build_nad83_harn.sh](https://raw.githubusercontent.com/OSGeo/proj-datumgrid/master/north-america/build_nad83_harn.sh) script.

* *hpgn.gsb (51 files)

### Continental USA VERTCON: NGVD (19)29 height to NAVD (19)88 height

*Source*: converted from [NOAA VERTCON .94 grids](https://www.ngs.noaa.gov/PC_PROD/VERTCON/)
*Converter*: vertcon_94_to_gtx.py
*Format*: GTX
*License*: Public Domain

Transform NGVD (19)29 height to NAVD (19)88 height. The .gtx files contain
the offset (in millimeter) to add to NGVD 29 heights to get NAVD 88 heights.

So with PROJ 5.2, to convert from NGVD 29 to NAVD 88, you need to a multiplication
by 0.001 to convert from millimeter to meter:
+proj=pipeline +step=vgridshift +multiplier=0.001 +grids=vertconw.gtx,vertconc.gtx,vertcone.gtx

The horizontal grid coordinates may be referenced to NAD27 or NAD83

* vertconc.gtx
* vertcone.gtx
* vertconw.gtx

### USA GEOID12B model: NAVD88/PRVD02/VIVD09/GUVD04/NMVD03/ASVD02 height to NAD83(2011)/NAD83(MA11)/NAD83(PA11)

*Source*: converted from [NOAA NGS .bin grids](https://www.ngs.noaa.gov/GEOID/GEOID12B/)
*Converter*: build_g2012.sh
*Format*: GTX
*License*: Public Domain

GEOID12B is a hybrid geoid height model that can transform between the relevant
vertical datum for each region to a NAD 83 ellipsoid height
See https://www.ngs.noaa.gov/GEOID/GEOID12B/faq_2012B.shtml for more details.

* g2012bu0.gtx : Conterminous US, referenced to NAD83(2011). Vertical datum: NAVD88
* g2012ba0.gtx : Alaska, referenced to NAD83(2011). Vertical datum: NAVD88
* g2012bh0.gtx : Hawaii, referenced to NAD83(PA11). Vertical datum: identical to that of USGG2012, which is offset by 50-60 cm from some Local Tidal bench mark values. 
* g2012bg0.gtx : Guam and Northern Mariana Islands, referenced to NAD83(PA11). Vertical datum: GUVD04/NMVD03
* g2012bs0.gtx : American Samoa, referenced to NAD83 (PA11). Vertical datum: ASVD02
* g2012bp0.gtx : Puerto Rico / U.S. Vigin Islands, referenced to NAD83(2011). Vertical datum: PRVD02 and VIVD09

## Included init-files

### Greenland: System definitions related to GR96

*Source*: [Agency for Data Supply and Efficiency](https://github.com/NordicGeodesy/NordicTransformations)  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Init files with various transformations and conversions of relevance in
Greenland. All definitions in the init file uses the local GR96 reference
frame as a pivot datum. More information can be found in the
[NordicTransformations](https://github.com/NordicGeodesy/NordicTransformations)
repository.

* GL
# PROJ-DATUMGRID-OCEANIA

The files in this package can be unpacked in the PROJ data directory.
For an installed PROJ this may be /usr/local/share/proj or /usr/share/proj
on unix style operating systems.

## About the Oceania PROJ resource package

The Oceania PROJ resource package is a collections of grids
and init-files that are not essential for the functionality of PROJ but
still of general interest to PROJ-users in the region specific to the
package.

## Included grids

### Australia: AGD66 -> GDA94

*Source*: [ICSM](https://github.com/icsm-au/transformation_grids)  
*Format*: NTv2  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Full national coverage two-way transformation file – see Appendix B
in the [GDA2020 Technical Manual](http://www.icsm.gov.au/gda/tech.html) for details.

* A66_National_13_09_01.gsb

### Australia: AGD66 -> GDA94

*Source*: [ICSM](https://github.com/icsm-au/transformation_grids)  
*Format*: NTv2  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Full coverage two-way transformation file in those states that
originally adopted AGD84 - West Australia, South Australia, Queensland
- see Appendix B in the [GDA2020 Technical Manual](http://www.icsm.gov.au/gda/tech.html)
for details.

* National_84_02_07_01.gsb

### Australia: GDA94 -> GDA2020

*Source*: [ICSM](https://github.com/icsm-au/transformation_grids)  
*Format*: NTv2  
*License*: [Creative Commons Attribution 4.0](https://creativecommons.org/licenses/by/4.0/)

Full national coverage two-way transformation files – See Sections 3.2
and 3.7 in the [GDA2020 Technical Manual](http://www.icsm.gov.au/gda/tech.html) for details.
Table 3.4 details which grid to use.

* GDA94_GDA2020_conformal.gsb
* GDA94_GDA2020_conformal_and_distortion.gsb

# PROJ-DATUMGRID-WORLD

The files in this package can be unpacked in the PROJ data directory.
For an installed PROJ this may be /usr/local/share/proj or /usr/share/proj
on unix style operating systems.

## About the World PROJ resource package

The World PROJ resource package is a collections of grids that are too large
to be in the general purpose proj-datumgrid package and not essential for
the functionality of PROJ, but still of general interest to PROJ-users.

## Included grids

### Vertical grid: EGM2008 geoid model

*Source*: [NGA](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html)  
*Format*: GTX
*License*: Public Domain

2.5 minute worldwide geoid undulation grid that transforms physical heights to WGS84 ellipsoidal heights.

This file has been produced by [GeographicLib](https://geographiclib.sourceforge.io/html/gravity.html)
using the EGM2008 gravity model and can be regenerated with the
[build_egm08_25_gtx.sh](https://raw.githubusercontent.com/OSGeo/proj-datumgrid/master/world/build_egm08_25_gtx.sh) script.

* egm08_25.gtx
