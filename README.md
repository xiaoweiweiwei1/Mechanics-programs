core codes:

double sigma3, sigma4,w=pi*pow(d,3)/32;
	sigma3 = sqrt(pow(F * L, 2)+m * m)/w;
	sigma4 = sqrt(pow(F * L, 2) + m * m*3/4)/w;
	if (sigma3 <= sigma)
		s = "safe";
	else s = "unsafe";
	if (sigma4 <= sigma)
		s1 = "safe";
	else s1 = "unsafe";
