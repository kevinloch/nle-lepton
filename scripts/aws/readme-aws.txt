###
### CAUTION: LARGE OR LONG RUNNING INSTANCES CAN BE VERY EXPENSIVE.   I AM NOT RESPONSIBLE FOR YOUR COSTS.
###
### Use responsibly and start with small instances/vcpus until you understand how it works
### and how much it will cost you.  Always manually set a termination date a few days in the future
### in case you forget they are running.   The default is one year!  Set up AWS billing alerts.
###

These files will help you run the lepton formula search program on as many AWS instances as you want.
The program is stateless and uses random startup seeding to be extremely parallelizable.  All results are
immediately uploaded to the server you configure in the source code so you can use cost-effective spot
instances without risk of data loss.

Instructions:

1. Start an instance of the same architecture type you intend to run in production (type 5c, 5d or 5n is recommended).
   Very little memory or disk space is needed (default AMI without block storage is fine) but highest performance cpu
   is desirable.

2. Install gcc with 'sudo bash; yum install gcc'.

3. Upload nle-lepton source code to server with scp or download from web server with curl, then compile with 'make'.

4. Copy the executable you compiled to your webserver (or another non-spot AWS instance running httpd) where start.sh will download it from.

5. Edit user-data.txt and start.sh with the desired url(s).

6. Edit nle-lepton.cfg with the desired arguments and put on your webserver.

7. Create a spot fleet request in the cheapest availability zone (currently US-EAST2 Ohio).   You can use a fixed instance type or specify number of vcpus.
   The default aws Linux AMI is fine.
   Make sure all of the instance types you let it choose from are the same architecture as what you compiled on (5* for example).   5c, 5d, 5n will all work
   with the same executable.  Copy and paste the contents of user-data.txt to the 'user data' field so it will start everything automatically.  Select 'maintain capacity'
   with 'terminate instances' option.   Set other options like security groups, expiration time as desired.

8. Once started the instances should grab start.sh and then grab nle-lepton.cfg and the executable from your webserver and start as many threads as there are vcpus on the server.

9. Watch your usage/costs carefully, and keep in mind there is a delay of about 12 hours before usage shows up in your aws billing dashboard.
   Even with the cheapest spot instances 100 vcpu's might cost around $.95/hour and that adds up fast ($23/day, $560/month)!

10. start.sh will check for Any change to nle-lepton.cfg once each minute.   If there are any changes it will kill all running threads, redownload the executable and restart
    the threads automatically.
