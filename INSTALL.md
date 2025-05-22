
```sh
git clone github.com/tornikeo/tblite
git clone http://github.com/tornikeo/tblite
code tblite/
```


```sh
pip install meson ninja ipykernel pandas numpy matplotlib seaborn
sudo apt-get install gfortran libblas-dev liblapack-dev -y
git config --global user.email tornikeonoprishvili@gmail.com
git config --global user.name tornikeo
meson setup build
meson compile -C build
meson test -C build/ hamiltonian --verbose
history 
```