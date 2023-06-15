source /home/jovyan/nemo/nemo_start.sh

magalie out=gal1mag_ic.snp ndisk=15000 nbulge=5000 nhalo=80000 mbulge=1/3 mhalo=16/3 rbulge=12 rhalo=12 
magalie out=gal2mag_ic.snp ndisk=15000 nbulge=5000 nhalo=80000 mbulge=1/3 mhalo=16/3 rbulge=12 rhalo=12 

snaprotate in=gal1mag_ic.snp out=gal1mag_rot.snp theta=60 spinvector=0.5,-0.8660254,0
snaprotate in=gal2mag_ic.snp out=gal2mag_rot.snp theta=60 spinvector=-0.5,0.8660254,0

snapshift in=gal1mag_rot.snp out=gal1mag_shift.snp rshift=[-4.6591016,5.65650751,0] vshift=[-0.16545881,-0.38691876,0]
snapshift in=gal2mag_rot.snp out=gal2mag_shift.snp rshift=[4.6591016,-5.65650751,0] vshift=[0.16545881,0.38691876,0]

snapstack in1=gal1mag_shift.snp in2=gal2mag_shift.snp out=col1mag.snp

gyrfalcON in=col1mag.snp out=col1mag_out.snp step=0.15 tstop=37.75 kmax=2.83 eps=0.015

snapscale in=col1mag_out.snp out=col1mag_final.snp mscale=0.1875 rscale=1/12 vscale=1.5