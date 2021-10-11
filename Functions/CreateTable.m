function CreateTable(d,v,a)
filename1='Acceleration-E(0%).txt';
filepath = 'F:\Plate to edit\text\';
file1=[filepath filename1];
fid1 = fopen(file1,'w');
fprintf(fid1,'%g\n' ,a);
fclose(fid1);
filename2='Displacement-E(0%).txt';
file2=[filepath filename2];
fid2 = fopen(file2,'w');
fprintf(fid2,'%g\n' ,d);
fclose(fid2);
filename3='Velocity-E(0%).txt';
file3=[filepath filename3];
fid3 = fopen(file3,'w');
fprintf(fid3,'%g\n' ,v);
fclose(fid3);
end