void decdeg2dms(float *deg,float *min,float *sec,float decdeg)
	{
	float decmin;
	*deg=(int) decdeg;
	decmin=(decdeg - *deg)*60.00;
	*min=(int) decmin;
	*sec=(decmin - *min)*60.00;
	}

