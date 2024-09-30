import os

# os.system( "rm splash.titles" )
os.system( "touch splash.titles" )

for i in range( 24 ):
    os.system( "echo \" snapshot %03d 600 kpc level 9 \" >> splash.titles" % i )