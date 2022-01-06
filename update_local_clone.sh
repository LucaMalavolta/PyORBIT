#!/bin/bash
git branch -m master legacy
git fetch origin
git branch -u origin/legacy legacy
git remote set-head origin -a
