Required env vars
```bash
export CC=/usr/bin/gcc-11
export CXX=/usr/bin/g++-11
export NVCC_PREPEND_FLAGS='-ccbin /usr/bin/g++-11'
export FC=/usr/bin/gfortran-13
```


```py
# Inputs

# NAO Number of spherical atomic orbitals in this basis set
# NSH number of shells in this basis set
# MAXL Maximum angular momentum of all basis functions, used to determine scratch size in integral calculation
mol = molecule struct
trans = zeros(10,10)
list = ... # adj list of atoms
bas = basis
h0 = ground state ham
selfenergy = zeros(NSH)... #diag hamilt elements zeros(N)
overlap = zeros(NAO,NAO) #, which elements interact

# Outputs
dpint = zeros(3, NAO, NAO) # dipole moment
qpint = zeros(6, NAO, NAO)
hamiltonian = zeros(NAO, NAO)

# Vars, inner
# int: 
#     iat, jat 
#     izp, jzp
#     itr, k, img, inl
#     ish, jsh,
#     is, jsm
#     ii, jj,
#     iao, jao, nao, ij
# float: 
#     rr, r2, vec(3), cutoff2, hij, shpoly, dtmpj, qtmpj
msao = [1, 3, 5, 7, 9, 11, 13] # This is a degree of symmetries?
stmp = zeros(msao[MAXL] ** 2)
dtmpi = zeros(3, msao[MAXL] ** 2)
qtmpi = zeros(6, msao[MAXL]**2)

## Others
#  h0.shpoly = zeros(all shells of all atoms, n unique atoms)
# Enhancement factor to scale the Hamiltonian elements ??
#  h0.hscale[jsh, ish, jzp, izp] 
#  h0.hscale[all shells, all shells, n unique atoms, n unique atoms] # How large is this in practice?
###  THE BIG LOOP

for iat in range(1, natoms)
    izp = # atom number
    is = # "shell space" offset for this atom
    inl = list.inl[iat] # neighbor map offset, which atom is adj to this
    for img in range(1, num_neighbors(iat)):
        jat = list.nlat[img + inl] # Index of neigh atom
        itr = list.nltr[img + inl] # "cell?" index of neigh atom
        jzp = id[jat] # atom num of neighbor
        js = bas.ish_at[jat] # SHell space offset for neighbor
        vec[:] = mol.xyz[:, iat] - mol.xyz[:, jat] - trans[:, itr] # Vector from us to neigh atom
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2 # Distance squared
        rr = sqrt(sqrt(r2) / (h0.rad[jzp] + h0.rad[izp])) # distance / sum radii of atoms, sqrt of that
        for ish in range(1, bas.nsh_id[izp]): # 1, to number of shells in this atom
            ii = bas.iao_sh[is + ish] # Index offset this shell in the atomic orbital space
            for jsh in range(1, bas.nsh_id[jzp]) # 1, to number of shells in neigh atom
                jj = bas.iao_sh[js + jsh] # get atom orbital space offset

                # TODO:
                multiple_cgto()

                # 1 + constant for [this shell and this atom number] * distance
                # x (1 + constant for [other shell, other atom number] * distance)
                shpoly = (1 + h0.shpoly[ish, izp] * rr) * (1 + h0.shpoly[jsh, jzp] * rr) # Just a scalar
                # Self energy for this shell + self enrgy for other shell
                hij = .5 * (selfenergy[is + ish] + selfenergy[js + jsh]) \
                         * h0.hscale[jsh, ish, jzp, izp] * shpoly # Just a scalar
                
                # get angular momentum jsh jzp basis function
                # 
                nao = msao[bas.cgto[jsh, jzp].ang ]
                for iao in range(1, msao[bas.cgto[ish,izp].ang]): # This atomic orb
                    for jao in range(1, nao): # Other atomic orb
                        ij = jao + nao*(iao - 1) # Indexing a matrix
                        # TODO
                        shift_operator()

                        # omp Atomic
                        overlap[jj+jao, ii+iao] += stmp[ij]

                        # omp atomic
                        for k in [1,2,3]:
                            dpint[k, jj+jao, ii+iao] += dtmpi[k,ij]
                        for k in range(1,6):
                            qpint[k,jj+jao,ii+iao] += qtmpi[k, ij]
                        hamiltonian[jj+jao,ii+iao] += stmp[ij] * hij

                        if iat != jat:
                            # If atoms are different, not self
                            overlap[ii+iao,jj+jao] += stmp[ij]
                            for k in [1,2,3]:
                                dpint[k, ii+iao, jj+jao] += dtmpj[k,ij]
                            for k in range(1,6):
                                qpint[k,ii+iao,jj+jao] += qtmpj[k]
                            hamiltonian[ii+iao,jj+jao] += stmp[ij] * hij

for iat in atoms: # For each atom 
    for ish in self_shells: # for all shels on that
        for jsh in self_shells: # To each shell
            multipole_cgto() # Some interaction integral
            for iao in msao(self): # All atomic orbitals in this shell?
                for jao in msao(self): # To all atomic orbitals...
                    ij = jao + nao * (iao - 1) # 
                    overlap += stmp[ij] # Overlap?
                    dpint[:, jj+jao,ii+iao] += dtmpi[:,ij]
                    qpint[:, jj+jao,ii+iao] += qtmpi[:,ij]
                    hamiltonian[jj+jao,ii+iao] += stmp[ij] * hij


```


```py
# multipole_cgto
# Contracted gaussian type basis function on Center i
# Contracted gaussian type basis function on Center j
# distance i to j
# r2 square distance
# vec(3) distance vector
# intcut  Maximum value of integral prefactor to consider
# overlap(msao(cgtoj%ang), msao(cgtoi%ang)) Overlap integrals for the given pair i  and j
# dpint(3, msao(cgtoj%ang), msao(cgtoi%ang)) Dipole moment integrals for the given pair i  and j
# qpint(6, msao(cgtoj%ang), msao(cgtoi%ang)) Quadrupole moment integrals for the given pair i  and j

# s3d(), d3d, q3d = zeros

for ip in range(1, nbasisi): # contraction len for this basis fun
    for jp in range(1, nbasisj):
        eab = cgtoi.alpha[ip] + cgtoj.alpha[jp]
        oab = 1/eab
        pre = exp(-est) * (sqrt(pi)**3) * sqrt(oab)**3
        cc = cgtoi.coeff[ip] * cgtoj.coeff[jp] * pre
        for mli in range(1, mlao[cgroj.ang]): # atomic orbitals in i
            for mlj in range(1, mlao[cgtoj.ang]): # atomic orbitals in j
                s3d[mlj,mli] += cc*val
                d3d[:,mlj,mli] += cc*dip
                q3d[:,mlj,mli] += cc*quad


```