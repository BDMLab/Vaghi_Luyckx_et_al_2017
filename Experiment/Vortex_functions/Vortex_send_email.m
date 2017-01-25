myaddress = 'data2server4storage@gmail.com';
mypassword = 'd4t4st0r4ge';

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

allfiles        = dir(fullfile(datafolder,'Vortex_ppt_*.mat'));
[blah,dx]       = sort([allfiles.datenum]);
relevantFile    = allfiles(dx(end)).name;

sendmail('data2server4storage@gmail.com', relevantFile, 'Vortex variance\n',{relevantFile});