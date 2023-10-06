echo "create new branch"
git init
git add CONFIG README.md *.sh script/** 
git commit -m "RelocaTE2 release"
git remote add origin https://github.com/JinfengChen/RelocaTE2_release.git
git push

echo "update"
git add CONFIG README.md *.sh script/**
git commit -m "update"
git remote add update https://github.com/JinfengChen/RelocaTE2_release.git
git push

echo "new rep for stajich lab relocate2"
git init
git pull https://github.com/stajichlab/RelocaTE2.git
git rm test_data.tar.gz
git commit -m "remove large test data"
git remote add update https://github.com/stajichlab/RelocaTE2.git
git push -u update master
