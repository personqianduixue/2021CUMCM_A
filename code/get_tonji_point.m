function single_tonji_point_xyz=get_tonji_point(xyz1,xyz2,xyz3)
n=8;
single_tonji_point_xyz=zeros(12,3);
AB_single_arrow=(xyz2-xyz1)./n;
AC_single_arrow=(xyz3-xyz1)./n;
D_xyz=xyz1+AB_single_arrow;
E_xyz=xyz1+AC_single_arrow;
DE_arrow=E_xyz-D_xyz;
k=1;
for i=1:4
    for j=0:(i-1)
       single_tonji_point_xyz(k,:)= xyz1+i*AB_single_arrow+j*DE_arrow;
       k=k+1;
    end
end
single_tonji_point_xyz(11,:)= xyz1+5*AB_single_arrow+2*DE_arrow;
single_tonji_point_xyz(12,:)= xyz1+5*AB_single_arrow+3*DE_arrow;

