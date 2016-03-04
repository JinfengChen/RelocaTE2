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

