#!/usr/bin/env nix-shell
#!nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/nixos-22.11.tar.gz
#! nix-shell -i bash -p ffmpeg
for i in ../data/motility/apng-videos/*.png; do ffmpeg -r 10 -i "$i" "${i%.*}.mkv"; done
