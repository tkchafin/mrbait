#!/bin/bash
make html
git add -A
git commit -m "update docs"
git push
