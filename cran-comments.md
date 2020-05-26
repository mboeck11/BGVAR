## Resubmission
This is a resubmission.
* R CMD check runs in 4mins. NO ERRORs, WARNINGs or NOTEs.
* Example now reduced to a toy example.
* Main reference added in the description. There are plenty of important references in the main estimation function as well as in the vignette.

## Resubmission
This is a resubmission. 

* fixed the issue with the invalid URL
* have added \donttest{} to most of the examples to reduce the time to run through
* the package involves quite demanding examples (to stick to examples from the literature) which need some minutes to compute, thus the vignette to build alone takes about 7-8mins
* checked via check_rhub()
	- Ubuntu Linux: 24mins
	- Debian Linux: 34mins
	- Windows Server: 38mins
	- Fedorea Linux: 41mins
	all with no ERRORs or WARNINGs and 2 NOTEs mentioned below.

## Test environments
* macOS Catalina 10.15.4, R 3.6.1

## R CMD check results
There wer no ERRORs or WARNINGs.
There are 2 NOTEs.

* checking CRAN incoming feasibility: possibly mis-spelled words in DESCRIPTION
- Autoregressions: plural of Autoregression
- NG: abbreviation of Normal-Gamma prior
- SSVS: abbreviation of Stochastic Search Variable Selection prior

Both abbreviations are explained in full detail, thus I do not assume this will cause problems.

* checking installed package size: installed size is 7.0Mb