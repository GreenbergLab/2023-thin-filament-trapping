(TeX-add-style-hook
 "scheme"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("mhchem" "version=4")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "mhchem"))
 :latex)

