program define epipop
	version 16.0
	`0'
end

program define simulate
	version 16.0
	gettoken model 0 : 0
	local cmd = "epipop_"+strlower(`"`model'"')
	`cmd' `0'
end

// end of file
