echo "====================="
echo "time series"
diff data/time_series.dat reference.out
echo "====================="
echo "powerhel_mag.dat"
diff data/powerhel_mag.dat reference_powerhel_mag.out
echo "====================="
echo "powerhel_kin.dat"
diff data/powerhel_kin.dat reference_powerhel_kin.out
echo "====================="
echo "power_mag.dat"
diff data/power_mag.dat reference_power_mag.out
echo "====================="
echo "power_kin.dat"
diff data/power_kin.dat reference_power_kin.out
echo "====================="
echo "power_krms.dat"
diff data/power_krms.dat reference_power_krms.out
echo "====================="
echo "power_lr.dat"
diff data/power_lr.dat reference_power_lr.out
echo "====================="
echo "poweru.dat"
diff data/poweru.dat reference_poweru.out
echo "====================="
echo "xyaverages.dat"
diff data/xyaverages.dat reference_xyaverages.out
echo "====================="
