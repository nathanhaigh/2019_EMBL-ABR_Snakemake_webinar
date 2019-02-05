This repository contains information for the EMBL-ABR webinar on "Snakemake: Repoducible and Scalable Bioifnormatic Workflows" presented by Nathan Watson-Haigh on 7th Feb 2019.

# Webinar Details

**Abstract:** Recent years have seen a groundswell of support in the bioscience community for improved reproducibility of data analyses. One solution used by bioinformaticians to achieve this is the use of the Snakemake workflow management system. It is a tool for creating reproducible and scalable data analyses.

This 1-hour webinar will take a look into Snakemake, cover the core concepts and examples required for making your first steps into reproducible Snakemake workflows. Building on this foundation, Nathan will introduce more advanced features which enable workflows to scale seamlessly to server, cluster, grid and cloud environments as well as different software execution environments.

**Date/Time:** Thursday 07 February 2019 [12:30-13:30 ACDT](https://www.timeanddate.com/worldclock/fixedtime.html?msg=EMBL-ABR+Snakemake+Webinar&iso=20190207T1230&p1=5&ah=1)

**Registration Page:** https://attendee.gotowebinar.com/register/7076869548376639747

**Presenter:** [Nathan Watson-Haigh](https://researchers.adelaide.edu.au/profile/nathan.watson-haigh)

**Video Link(s):** TBA following the webinar

# Data for the Webinar

For the purpose of demonstrating a Snakemake workflow in reasonable time, we will be working on a subset of public WGS data from wheat. Specifically, we will be looking at data from a small (58Kbp) region on the long arm of chromosome 4A (chr4A_part2:235500000-235558000). This region contains 2 genes affected by a deletion which is present in some wheat accessions:

![alt text](img/chr4A_Wx-B1_Null_region.png  "chr4A_part2:235500000-235558000")

DAWN URL: [chr4A_part2:235500000-235558000](http://crobiad.agwine.adelaide.edu.au/dawn/jbrowse/?loc=chr4A_part2%3A235500000..235558000&tracks=IWGSC_v1.0_HC_genes%2CAlsen_snpcoverage%2CRAC875_snpcoverage%2CYitpi_snpcoverage%2CPastor_snpcoverage%2CWyalkatchem_snpcoverage%2CWestonia_snpcoverage%2CACBarrie_snpcoverage%2CVolcanii_snpcoverage%2CBaxter_snpcoverage%2CChara_snpcoverage%2CDrysdale_snpcoverage%2CH45_snpcoverage%2CXiaoyan_snpcoverage%2CKukri_snpcoverage%2CGladius_snpcoverage%2CExcalibur_snpcoverage)
