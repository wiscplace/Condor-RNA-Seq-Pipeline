#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: sendEmail.py

@Purpose: Send email with message
          
@Input  :  Text Message            
                     
@Output : email

UW DOIT :   smtp.wiscmail.wisc.edu
SMTP Server: relay.mail.wisc.edu or smtp.wiscmail.wisc.edu
SMTP Port: 25

@author: Mike Place
@Date:   3/18/2016
"""
import os
import smtplib
import sys

def send( message ):
    smtpObj = smtplib.SMTP('smtp.wiscmail.wisc.edu',25)
    text = 'Subject: Condor Job status\n\n' + message
    smtpObj.sendmail('mplace@wisc.edu', 'mplace@wisc.edu', text )
    smtpObj.quit()    

#def main():
#    """
#    Main 
#    """
#    message = "RNA-Seq processing complete, have a jolly day" # text


if __name__ == "__main__":
    main()
 
