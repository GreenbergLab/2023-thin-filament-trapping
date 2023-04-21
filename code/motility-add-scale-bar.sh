#!/usr/bin/env nix-shell
#!nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/nixos-22.11.tar.gz
#! nix-shell -i bash -p ffmpeg
ffmpeg -i ../data/motility/apng-videos/to-crop/biotin-actin-pCa-4-labeled-cropped.mkv \
-vf "drawbox=x=iw-w-10:y=ih-h-10:w=95:h=5:color=white,\
drawtext=font='Helvetica':text='10 microns':fontcolor=white:fontsize=20:box=1:boxcolor=black:boxborderw=0:x=w-110:y=h-40" \
../data/motility/apng-videos/to-crop/biotin-actin-pCa-4-labeled-cropped-scale-bar.mkv
