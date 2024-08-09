ps -ef | grep visualServer | grep -v grep| awk '{print "kill -9 " $2}'| sh
ps -ef | grep python | grep -v grep| awk '{print "kill -9 " $2}'| sh

