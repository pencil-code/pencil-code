amin=-0.1
amax=+6.5
pc_read_grid,obj=grid,/trim
rvid_line,'phi',proc=0,min=amin,max=amax,map=map2,tt=tt,xaxis=x
contour,transpose(map2),nlev=30,/fil,tt,grid.x,xtit='!8t!6',ytit='!8x!6'
END
