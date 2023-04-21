#!/usr/bin/env nix-shell
#!nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/nixos-22.11.tar.gz
#! nix-shell -i bash -p ffmpeg
for i in ../data/motility/apng-videos/to-crop/*-labeled.mkv; do ffmpeg -i "$i" -filter:v "crop=in_w*0.66:in_h*0.66:0:0" "${i%.*}"-cropped.mkv; done
