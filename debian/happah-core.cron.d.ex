#
# Regular cron jobs for the happah-core package
#
0 4	* * *	root	[ -x /usr/bin/happah-core_maintenance ] && /usr/bin/happah-core_maintenance
