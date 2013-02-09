{
    NX=$1;
    NY=$2;
    NZ=$3;
#    print NX,NY,NZ;

    for(k=0;k<NZ;k++)
    for(j=0;j<NY;j++)
    for(i=0;i<NX;i++){
	id=i+NX*j+NX*NY*k;
	ip=i+1;
	im=i-1;
        jp=j+1;
        jm=j-1;
        kp=k+1;
        km=k-1;
	if(ip<NX && ip>=0){id_ip=ip+NX*j+NX*NY*k;}else{id_ip="-1";}
	if(im<NX && im>=0){id_im=im+NX*j+NX*NY*k;}else{id_im="-1";}
	if(jp<NY && jp>=0){id_jp=i+NX*jp+NX*NY*k;}else{id_jp="-1";}
	if(jm<NY && jm>=0){id_jm=i+NX*jm+NX*NY*k;}else{id_jm="-1";}
	if(kp<NZ && kp>=0){id_kp=i+NX*j+NX*NY*kp;}else{id_kp="-1";}
	if(km<NZ && km>=0){id_km=i+NX*j+NX*NY*km;}else{id_km="-1";}       

#	print id,i,j,k,ip,id_ip,id_im, id_jp, id_jm, id_kp,id_km;
	printf "%e %e %e %d %d %d %d %d %d\n", i,j,k,id_ip,id_im, id_jp, id_jm, id_kp,id_km; 
    }

}
