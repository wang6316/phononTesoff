%main routing calculating phonon structrue of the superlattice using
%Tersoff potential

clear
clc

tic

%build a framework for the superlattice
n1=1;%number of 8-atom cells in well region
n2=1;%number of 8-atom cells in barrier region
ndl=n1+n2;%number of double layers
na=ndl*8;%number of atoms in each supercell
nb=3*na;%number of branches of phonon for each wavevector

%lattice constant of the supercell, scale on "a"
xlatt=1;
ylatt=ndl;
zlatt=1;

%mass in atomic mass unit "g/mol"
mAl=26.98154;
mIn=114.82;
mAs=74.9216;
mSb=121.75;


%look up table for atoms in superlattice
mass=zeros(na,1);
species=zeros(na,1);
%assign atoms to each site
for n=1:n1
    mass(1+(n-1)*8)=mIn;
    mass(3+(n-1)*8)=mIn;
    mass(5+(n-1)*8)=mIn;
    mass(7+(n-1)*8)=mIn;
    
    species(7+(n-1)*8)=1;
    species(1+(n-1)*8)=1;
    species(3+(n-1)*8)=1;
    species(5+(n-1)*8)=1;
    
    mass(2+(n-1)*8)=mAs;
    mass(4+(n-1)*8)=mAs;
    mass(6+(n-1)*8)=mAs;
    mass(8+(n-1)*8)=mAs;
    
    species(2+(n-1)*8)=2;
    species(4+(n-1)*8)=2;
    species(6+(n-1)*8)=2;
    species(8+(n-1)*8)=2;
end

for n=(n1+1):(n1+n2)
    mass(1+(n-1)*8)=mAl;
    mass(3+(n-1)*8)=mAl;
    mass(5+(n-1)*8)=mAl;
    mass(7+(n-1)*8)=mAl;
    
    species(1+(n-1)*8)=3;
    species(3+(n-1)*8)=3;
    species(5+(n-1)*8)=3;
    species(7+(n-1)*8)=3;    
    
    mass(2+(n-1)*8)=mSb;
    mass(4+(n-1)*8)=mSb;
    mass(6+(n-1)*8)=mSb;
    mass(8+(n-1)*8)=mSb;
    
    species(1+(n-1)*8)=4;
    species(3+(n-1)*8)=4;
    species(5+(n-1)*8)=4;
    species(7+(n-1)*8)=4;
end

%mass matrix to normalize dynamical matrix
mass3dn=zeros(3*na);
for i=1:na
    mass3dn(3*i-2,3*i-2)=1/sqrt(mass(i));
    mass3dn(3*i-1,3*i-1)=1/sqrt(mass(i));
    mass3dn(3*i,3*i)=1/sqrt(mass(i));
end

%relative coordinates of atoms in a 8-atom cell
%(default values, we can include lattice distortion here)
%(length scale on "a" of ffc)
tao=cell(8,1);
tao{1}=[0,0,0];
tao{2}=tao{1}+[1/4,1/4,1/4];
tao{3}=[1/2,1/2,0];
tao{4}=tao{3}+[1/4,1/4,1/4];
tao{5}=[1/2,0,1/2];
tao{6}=tao{5}+[1/4,1/4,1/4];
tao{7}=[0,1/2,1/2];
tao{8}=tao{7}+[1/4,1/4,1/4];

%coodinates in supercell
R=cell(na,1);
for i=1:ndl
    for p=1:8
        nn=8*(i-1)+p;
        R{nn}=[0,(i-1),0]+tao{p};
    end
end

%evaluate interatomic force constant using Tersoff potential 
%fill in informations about 4 nearest neighbors and 12 next nearest neighbors 

%assign lattice constant(in angstrom)
a=5.620;

nn=cell(na,1);
nntype=cell(na,1);
nnlatt=cell(na,1);%lattice index of nearest neighbor

for i=1:na%loop i runs over all atoms of the supercell 
    itype=species(i);
    imass=mass(i);
    nn{i}=cell(16,1);
    nnlatt{i}=cell(16,1);
    nntype{i}=zeros(16,1);
    inn=0;
    for ix=-1:1:1
        for iy=-1:1:1
            for iz=-1:1:1               
                for j=1:na
                    tempR=R{j}+[xlatt*ix,ylatt*iy,zlatt*iz];
                    distance=norm(tempR-R{i});                   
                    if (distance>0.4 && distance<0.75)
                        inn=inn+1;
                        nn{i}{inn}=tempR;
                        nnlatt{i}{inn}=[ix,iy,iz];
                        nntype{i}(inn)=species(j);
                    end
                end
            end
        end
    end
end
% 
% %finding force constant using Tersoff potential
% for i=na%loop for atoms in a unit cell
%     icoord=R{i};
%     for j=1:4 %loops for 4 nearest neighbors of atom i, within cutoff for current parameters
%         jcoord=nn{}
        
        
        
        
        
        
                
            


    



    






    
