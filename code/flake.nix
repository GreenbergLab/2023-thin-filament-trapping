{
  description = "Developmental environment for publication from the Greenberg Lab: 2023-thin-filament-trapping";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      #nixpkgs
      pkgs = nixpkgs.legacyPackages.${system};

    tex = (pkgs.texlive.combine {
      inherit (pkgs.texlive) scheme-basic mhchem;
      });

    in {
      devShells.default = pkgs.mkShell {
        nativeBuildInputs = [ pkgs.bashInteractive ];
        buildInputs = with pkgs; [ R
                                   tex
                                   ghostscript
                                   rPackages.ggplot2
                                   rPackages.purrr
                                   rPackages.readr
                                   rPackages.tidyr
                                   rPackages.dplyr
                                   rPackages.ggtext
                                   rPackages.here
                                   rPackages.readxl
                                   rPackages.data_table
                                   rPackages.ggpubr
                                   rPackages.cowplot
                                   rPackages.magick
                                   rPackages.rsvg
                                   rPackages.pdftools];
      };
    });
}
