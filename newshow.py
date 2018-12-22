#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
plt.switch_backend('agg');
path=sys.argv[1];
raw_fraction="0 0.05 0.07 0.08 0.085 0.09 0.1 0.11 0.12 0.13 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.8 1.0";
fraction=[];
for i in raw_fraction.split():
	fraction.append(float(i));
l_x=[];
l_y=[];
l_z=[];
d_x_zr=[];
d_y_zr=[];
d_z_zr=[];
d_x_ba=[];
d_y_ba=[];
d_z_ba=[];
d_x_ca=[];
d_y_ca=[];
d_z_ca=[];
d_scalar_zr=[];
d_scalar_ba=[];
d_scalar_ca=[];
tilt_angle=[];
px=[];
py=[];
pz=[];
pall=[];
dic=[];
epsilonx=[];
epsilony=[];
epsilonz=[];
print fraction
for i in fraction:
	strtemp=path+"/BaCa"+str(format(i,".2f"))+"/result.txt";
	f=open(strtemp,"r");
	la_const=f.readline();
	la_const=f.readline();
	la_const=la_const.split();
	temp=[float(la_const[0]),float(la_const[1]),float(la_const[2])];
	temp.sort();
	l_x.append(temp[0]);
	l_y.append(temp[1]);
	l_z.append(temp[2]);	
	volume=temp[0]*temp[1]*temp[2];
	da_const=f.readline();
	da_const=f.readline();
	da_const=da_const.split();
	temp=[float(da_const[0]),float(da_const[1]),float(da_const[2])];
	temp.sort();
	d_x_zr.append(temp[0]);
	d_y_zr.append(temp[1]);
	d_z_zr.append(temp[2]);
	da_const=f.readline();
	da_const=f.readline();
	da_const=da_const.split();
	temp=[float(da_const[0]),float(da_const[1]),float(da_const[2])];
	temp.sort();
	d_x_ba.append(temp[0]);
	d_y_ba.append(temp[1]);
	d_z_ba.append(temp[2]);
	da_const=f.readline();
	da_const=f.readline();
	d_x_ca.append(temp[0]);
	d_y_ca.append(temp[1]);
	d_z_ca.append(temp[2]);
	da_const=f.readline();
	da_const=f.readline();
	da_const=float(da_const)
	tilt_angle.append(da_const);
	da_const=f.readline();
	da_const=f.readline();
	d_scalar_zr.append(float(da_const));
	da_const=f.readline();
	da_const=f.readline();
	d_scalar_ba.append(float(da_const));
	da_const=f.readline();
	da_const=f.readline();
	d_scalar_ca.append(float(da_const));
	da_const=f.readline();
	da_const=f.readline();
	da_const=da_const.split();
	tempall=[float(da_const[0]),float(da_const[1]),float(da_const[2])];
	da_const=f.readline();
	da_const=f.readline();
	da_const=da_const.split();
	temp=[float(da_const[0]),float(da_const[1]),float(da_const[2])];
	dic={tempall[i]:temp[i] for i in range(len(temp))};
	rank=sorted(dic);
	epsilonx.append(dic[rank[0]]*volume);
	epsilony.append(dic[rank[1]]*volume);
	epsilonz.append(dic[rank[2]]*volume);
l_v=np.multiply(l_x,l_y);
l_v=np.multiply(l_v,l_z);
l_v=[i*8 for i in l_v];
fig=plt.figure(1);
plt.plot(fraction,l_v,"ro--",label="c")
plt.legend();
plt.xlabel("different fraction of Calcium");
plt.ylabel("crystal volume 2*2*2");
plt.savefig("crystal_volume.png");
fig=plt.figure(2);
d_all_zr=np.sqrt(np.square(d_x_zr)+np.square(d_y_zr)+np.square(d_z_zr));
plt.plot(fraction,d_z_zr,"ro--",fraction,d_x_zr,"bs--",fraction,d_y_zr,"g^--",fraction,d_all_zr,"+--");
plt.legend(["dz","dx","dy","dall"]);
plt.xlabel("different fraction of calcium");
plt.ylabel("displacement of Zr");
plt.savefig("displacement_zr");
fig=plt.figure(3);
plt.plot(fraction,l_x,"ro--",fraction,l_y,"bs--",fraction,l_z,"g^--");
plt.legend(["lattice_x","lattice_y","lattice_z"]);
plt.xlabel("different fraction of calcium");
plt.ylabel("lattice constant of the crystal");
plt.savefig("lattice_const");
fig=plt.figure(4);
d_all_ba=np.sqrt(np.square(d_x_ba)+np.square(d_y_ba)+np.square(d_z_ba));
plt.plot(fraction,d_z_ba,"ro--",fraction,d_x_ba,"bs--",fraction,d_y_ba,"g^--",fraction,d_all_ba,"+--");
plt.legend(["dz","dx","dy","dall"]);
plt.xlabel("different fraction of calcium");
plt.ylabel("displacement of ba");
plt.savefig("displacement_ba");
fig=plt.figure(5);
d_all_ca=np.sqrt(np.square(d_x_ca)+np.square(d_y_ca)+np.square(d_z_ca));
plt.plot(fraction,d_z_ca,"ro--",fraction,d_x_ca,"bs--",fraction,d_y_ca,"g^--",fraction,d_all_ca,"+--");
plt.legend(["dz","dx","dy","dall"]);
plt.xlabel("different fraction of calcium");
plt.ylabel("displacement of ca");
plt.savefig("displacement_ca");
fig=plt.figure(6);
d_all_ca=np.sqrt(np.square(d_x_ca)+np.square(d_y_ca)+np.square(d_z_ca));
plt.plot(fraction,d_all_ca,"ro--",fraction,d_all_ba,"bs--",fraction,d_all_zr,"g^--");
plt.legend(["Ca","Ba","Zr"]);
plt.xlabel("different fraction of calcium");
plt.ylabel("displacement of all the atoms");
plt.savefig("displacement_all");
fig=plt.figure(7);
plt.plot(fraction,tilt_angle);
plt.legend("average tilt angles");
plt.xlabel("different fraction of calcium");
plt.ylabel("tilt angle of O6 cage");
plt.savefig("tilt-angle");
final=open("final.txt","w")
fig=plt.figure(8);
plt.plot(fraction,d_scalar_zr,"ro--",fraction,d_scalar_ba,"bs--",fraction,d_scalar_ca,"g^--");
plt.legend(["Zr displacement","ba displacement","Ca displacement"]);
plt.xlabel("calcium concentration");
plt.ylabel("displacement/A");
plt.savefig("displacement_scalar.png");
fig=plt.figure(9);
print epsilonx
plt.plot(fraction,epsilonx,"ro--",fraction,epsilony,"bs--",fraction,epsilonz,"g^--");
plt.legend(["epsilon_x","epsilon_y","epsilon_z"]);
plt.xticks(np.arange(0,1,0.1))
plt.xlabel("concentration");
plt.ylabel("epsilon");
plt.savefig("epsilon.png")