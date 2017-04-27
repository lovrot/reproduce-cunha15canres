#!/bin/bash

git branch -D gh-pages  # delete old branch -- no reason to version control gh-pages
git push origin :gh-pages

cd /tmp/
git clone ~/git/reproduce-cunha15canres

cd /tmp/reproduce-cunha15canres/
git checkout --orphan gh-pages  # create new
git rm -rf .

cp ~/git/reproduce-cunha15canres/reports/reproduce-cunha15canres.html ./index.html
cp ~/git/reproduce-cunha15canres/reports/addendum-train-illustrative-models.html ./

git add index.html
git commit -m "Update gh-page"
git push origin gh-pages

cd ..
rm -fR reproduce-cunha15canres

cd ~/git/reproduce-cunha15canres

git push origin gh-pages
