"""
Code to go into my gmail account and get the Canadian Fraser River data
they send me every day.  It only saves the attachments that have not already
been saved (nice!).

We have emails beginning 2014.02.12.

Performance: about a half second per file, or 5 min for 600 files,
but it only takes a couple of seconds in practice becase it just looks at
mails that we haven't already pulled the attachement from (using ndays).

Adapted from:
http://stackoverflow.com/questions/7596789/
    downloading-mms-emails-sent-to-gmail-using-python
"""

import os; import sys; pth = os.path.abspath('../../alpha')
if pth not in sys.path: sys.path.append(pth)
import Lfun; reload(Lfun)
Ldir = Lfun.Lstart()

import email, imaplib, os
from datetime import datetime

# get account info
g_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/gmail_pm_2015.10.17.csv')

outdir = Ldir['data'] + 'rivers/data_fraser_raw/'
# directory in which to save attachments

# what is the last day we have?
outlist = os.listdir(outdir)
dt_last = datetime.strptime(outlist[-1],'Parker%Y%m%d.csv')
dt_now = datetime.now()
how_long_ago = dt_now - dt_last
ndays = how_long_ago.days
# use this to speed up the search for new attachments

# connecting to the gmail imap server
m = imaplib.IMAP4_SSL("imap.gmail.com")
m.login(g_dict['username'], g_dict['password'])
m.select('LiveOcean/Rivers_Canada')
# here you a can choose a mail box like INBOX instead
# use m.list() to get all the mailboxes

# you could filter using the IMAP rules here
# (check http://www.example-code.com/csharp/imap-search-critera.asp)
resp, items = m.search(None, "ALL")

items = items[0].split()
# getting the mails id

for emailid in items[-(ndays + 3):]:
    resp, data = m.fetch(emailid, "(RFC822)")
    # fetching the mail, "`(RFC822)`" means "get the whole stuff",
    # but you can ask for headers only, etc.
    
    email_body = data[0][1]
    # getting the mail content
    
    mail = email.message_from_string(email_body)
    # parsing the mail content to get a mail object
    #
    # Note: you could get the datetime of the message using
    # from datetime import datetime
    # dt = datetime.strptime(mail['Date'][:-6],'%d %b %Y %H:%M:%S')
    # You have to drop the last 6 characters becasue they are the timezone
    # offset (e.g. -0800) but parsing this is not part of python 2.7.

    #Check if any attachments at all
    if mail.get_content_maintype() != 'multipart':
        continue

    # We use walk to create a generator so we can iterate on the parts
    # and forget about the recursive headache.
    for part in mail.walk():
        
        # multipart are just containers, so we skip them
        if part.get_content_maintype() == 'multipart':
            continue

        # is this part an attachment ?
        if part.get('Content-Disposition') is None:
            continue

        filename = part.get_filename()
        counter = 1

        # if there is no filename, we create one with a counter to avoid duplicates
        if not filename:
            filename = 'part-%03d%s' % (counter, 'bin')
            counter += 1

        att_path = os.path.join(outdir, filename)

        #Check if its already there
        if not os.path.isfile(att_path):
            # finally write the stuff
            print("[" + mail["From"] + "] :" + mail["Subject"] + ' : ' + filename)
            fp = open(att_path, 'wb')
            fp.write(part.get_payload(decode=True))
            fp.close()