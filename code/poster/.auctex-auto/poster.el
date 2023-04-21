(TeX-add-style-hook
 "poster"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("beamer" "final")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("beamerposter" "size=custom" "width=127" "height=106.68" "scale=1.15") ("caption" "labelfont=bf") ("natbib" "numbers") ("microtype" "patch=none")))
   (add-to-list 'LaTeX-verbatim-environments-local "semiverbatim")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "href")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "beamer"
    "beamer10"
    "fontenc"
    "lmodern"
    "beamerposter"
    "graphicx"
    "booktabs"
    "doi"
    "caption"
    "natbib"
    "microtype"
    "tikz"
    "pgfplots"
    "anyfontsize"
    "textgreek")
   (TeX-add-symbols
    "separatorcolumn"
    "translate")
   (LaTeX-add-lengths
    "sepwidth"
    "colwidth"))
 :latex)

