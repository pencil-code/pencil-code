Kernel impose_floors_and_ceilings()
{
	if(AC_density_floor__mod__density > 0.)
	{
		if(AC_enum_density_floor_profile__mod__density == enum_uniform_string)
		{
			LNRHO[vertexIdx.x][vertexIdx.y][vertexIdx.z] = max(AC_density_floor_log__mod__density,LNRHO[vertexIdx.x][vertexIdx.y][vertexIdx.z])
		}
	}
}
