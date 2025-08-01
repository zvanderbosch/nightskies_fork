#!/usr/bin/env python3

# This file is part of the Astrometry.net suite.
# Copyright 2009 Dustin Lang
# Licensed under a 3-clause BSD style license - see LICENSE
# https://github.com/dstndstn/astrometry.net

from __future__ import print_function
import os
import sys
import time
import json

# py3
from urllib.parse import urlencode, quote
from urllib.request import urlopen, Request

# local imports
import printcolors as pc

# Define print status prefix
scriptName = 'client.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '


def json2python(data):
    try:
        return json.loads(data)
    except:
        pass
    return None
python2json = json.dumps

class MalformedResponse(Exception):
    pass
class RequestError(Exception):
    pass

class Client(object):
    default_url = 'https://nova.astrometry.net/api/'

    def __init__(self,
                 apiurl = default_url):
        self.session = None
        self.apiurl = apiurl

    def get_url(self, service):
        return self.apiurl + service

    def send_request(self, service, args={}, file_args=None):
        '''
        service: string
        args: dict
        '''
        if self.session is not None:
            args.update({ 'session' : self.session })
        json = python2json(args)
        url = self.get_url(service)

        # If we're sending a file, format a multipart/form-data
        if file_args is not None:
            import random
            boundary_key = ''.join([random.choice('0123456789') for i in range(19)])
            boundary = '===============%s==' % boundary_key
            headers = {'Content-Type':
                       'multipart/form-data; boundary="%s"' % boundary}
            data_pre = (
                '--' + boundary + '\n' +
                'Content-Type: text/plain\r\n' +
                'MIME-Version: 1.0\r\n' +
                'Content-disposition: form-data; name="request-json"\r\n' +
                '\r\n' +
                json + '\n' +
                '--' + boundary + '\n' +
                'Content-Type: application/octet-stream\r\n' +
                'MIME-Version: 1.0\r\n' +
                'Content-disposition: form-data; name="file"; filename="%s"' % file_args[0] +
                '\r\n' + '\r\n')
            data_post = (
                '\n' + '--' + boundary + '--\n')
            data = data_pre.encode() + file_args[1] + data_post.encode()

        else:
            # Else send x-www-form-encoded
            data = {'request-json': json}
            data = urlencode(data)
            data = data.encode('utf-8')
            headers = {}

        # Generate the request
        request = Request(url=url, headers=headers, data=data)
        
        # Try sending the request
        delay = 1.0      # Seconds between retry attempts
        max_retries = 20 # Max number of retries
        for attempts in range(max_retries):
            try:
                f = urlopen(request)
                txt = f.read()
                result = json2python(txt)
                stat = result.get('status')
                if stat == 'error':
                    errstr = result.get('errormessage', '(none)')
                    raise RequestError('server error message: ' + errstr)
                return result
            except:
                if attempts < max_retries-1:
                    time.sleep(delay)
                    delay += 0.5
                else:
                    print(
                        f'{PREFIX}(Status = {pc.RED}FAILED{pc.END}  ) '
                        f'Exceeded max retries while sending {service} request'
                    )
                    sys.exit(-1)


    def login(self, apikey):
        args = { 'apikey' : apikey }
        result = self.send_request('login', args)
        sess = result.get('session')
        if not sess:
            raise RequestError('no session in result')
        self.session = sess

    def _get_upload_args(self, **kwargs):
        args = {}
        for key,default,typ in [('allow_commercial_use', 'd', str),
                                ('allow_modifications', 'd', str),
                                ('publicly_visible', 'y', str),
                                ('scale_units', None, str),
                                ('scale_type', None, str),
                                ('scale_lower', None, float),
                                ('scale_upper', None, float),
                                ('scale_est', None, float),
                                ('scale_err', None, float),
                                ('center_ra', None, float),
                                ('center_dec', None, float),
                                ('parity',None,int),
                                ('radius', None, float),
                                ('downsample_factor', None, int),
                                ('positional_error', None, float),
                                ('tweak_order', None, int),
                                ('crpix_center', None, bool),
                                ('invert', None, bool),
                                ('image_width', None, int),
                                ('image_height', None, int),
                                ('x', None, list),
                                ('y', None, list),
                                ('album', None, str),
                                ]:
            if key in kwargs:
                val = kwargs.pop(key)
                val = typ(val)
                args.update({key: val})
            elif default is not None:
                args.update({key: default})
        #print('Upload args:', args)
        return args

    def url_upload(self, url, **kwargs):
        args = dict(url=url)
        args.update(self._get_upload_args(**kwargs))
        result = self.send_request('url_upload', args)
        return result

    def upload(self, fn=None, **kwargs):
        args = self._get_upload_args(**kwargs)
        file_args = None
        if fn is not None:
            try:
                f = open(fn, 'rb')
                file_args = (fn, f.read())
            except IOError:
                print(f'{PREFIX}File {fn} does not exist')
                raise
        return self.send_request('upload', args, file_args)

    def submission_images(self, subid):
        result = self.send_request('submission_images', {'subid':subid})
        return result.get('image_ids')


    def myjobs(self):
        result = self.send_request('myjobs/')
        return result['jobs']

    def job_status(self, job_id, justdict=False):
        result = self.send_request('jobs/%s' % job_id)
        if justdict:
            return result
        stat = result.get('status')
        if stat == 'success':
            result = self.send_request('jobs/%s/calibration' % job_id)
            print(f'{PREFIX}Calibration:', result)
            result = self.send_request('jobs/%s/tags' % job_id)
            print(f'{PREFIX}Tags:', result)
            result = self.send_request('jobs/%s/machine_tags' % job_id)
            print(f'{PREFIX}Machine Tags:', result)
            result = self.send_request('jobs/%s/objects_in_field' % job_id)
            print(f'{PREFIX}Objects in field:', result)
            result = self.send_request('jobs/%s/annotations' % job_id)
            print(f'{PREFIX}Annotations:', result)
            result = self.send_request('jobs/%s/info' % job_id)
            print(f'{PREFIX}Calibration:', result)

        return stat

    def annotate_data(self,job_id):
        """
        :param job_id: id of job
        :return: return data for annotations
        """
        result = self.send_request('jobs/%s/annotations' % job_id)
        return result
		
    def calibrate_data(self,job_id):
        """
        :param job_id: id of job
        :return: parity,orientation,pixscale,radius,ra, and dec
        """
        result = self.send_request('jobs/%s/calibration' % job_id)
        return result

    def sub_status(self, sub_id, justdict=False):
        result = self.send_request('submissions/%s' % sub_id)
        if justdict:
            return result
        return result.get('status')

    def jobs_by_tag(self, tag, exact):
        exact_option = 'exact=yes' if exact else ''
        result = self.send_request(
            'jobs_by_tag?query=%s&%s' % (quote(tag.strip()), exact_option),
            {},
        )
        return result

if __name__ == '__main__':
    #print("Running with args %s"%sys.argv)
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--server', dest='server', default=Client.default_url,
                      help='Set server base URL (eg, %default)')
    parser.add_option('--apikey', '-k', dest='apikey',
                      help='API key for Astrometry.net web service; if not given will check AN_API_KEY environment variable')
    parser.add_option('--upload', '-u', dest='upload', help='Upload a file')
    parser.add_option('--wait', '-w', dest='wait', action='store_true', help='After submitting, monitor job status')
    parser.add_option('--wcs', dest='wcs', help='Download resulting wcs.fits file, saving to given filename; implies --wait if --urlupload or --upload')
    parser.add_option('--newfits', dest='newfits', help='Download resulting new-image.fits file, saving to given filename; implies --wait if --urlupload or --upload')
    parser.add_option('--corr', dest='corr', help='Download resulting corr.fits file, saving to given filename; implies --wait if --urlupload or --upload')
    parser.add_option('--kmz', dest='kmz', help='Download resulting kmz file, saving to given filename; implies --wait if --urlupload or --upload')
    parser.add_option('--annotate','-a',dest='annotate',help='store information about annotations in give file, JSON format; implies --wait if --urlupload or --upload')
    parser.add_option('--calibrate','-C',dest='calibrate',help='store information about calibration in given file, JSON format; implies --wait if --urlupload or --upload')
    parser.add_option('--urlupload', '-U', dest='upload_url', help='Upload a file at specified url')
    parser.add_option('--scale-units', dest='scale_units',
                      choices=('arcsecperpix', 'arcminwidth', 'degwidth', 'focalmm'), help='Units for scale estimate')
    #parser.add_option('--scale-type', dest='scale_type',
    #                  choices=('ul', 'ev'), help='Scale bounds: lower/upper or estimate/error')
    parser.add_option('--scale-lower', dest='scale_lower', type=float, help='Scale lower-bound')
    parser.add_option('--scale-upper', dest='scale_upper', type=float, help='Scale upper-bound')
    parser.add_option('--scale-est', dest='scale_est', type=float, help='Scale estimate')
    parser.add_option('--scale-err', dest='scale_err', type=float, help='Scale estimate error (in PERCENT), eg "10" if you estimate can be off by 10%')
    parser.add_option('--ra', dest='center_ra', type=float, help='RA center')
    parser.add_option('--dec', dest='center_dec', type=float, help='Dec center')
    parser.add_option('--radius', dest='radius', type=float, help='Search radius around RA,Dec center')
    parser.add_option('--downsample', dest='downsample_factor', type=int, help='Downsample image by this factor')
    parser.add_option('--positional_error', dest='positional_error', type=float, help='How many pixels a star may be from where it should be.')
    parser.add_option('--parity', dest='parity', choices=('0','1'), help='Parity (flip) of image')
    parser.add_option('--tweak-order', dest='tweak_order', type=int, help='SIP distortion order (default: 2)')
    parser.add_option('--crpix-center', dest='crpix_center', action='store_true', default=None, help='Set reference point to center of image?')
    parser.add_option('--invert', action='store_true', default=None, help='Invert image before detecting sources -- for white-sky, black-stars images')
    parser.add_option('--image-width', type=int, help='Set image width for x,y lists')
    parser.add_option('--image-height', type=int, help='Set image height for x,y lists')
    parser.add_option('--album', type=str, help='Add image to album with given title string')
    parser.add_option('--sdss', dest='sdss_wcs', nargs=2, help='Plot SDSS image for the given WCS file; write plot to given PNG filename')
    parser.add_option('--galex', dest='galex_wcs', nargs=2, help='Plot GALEX image for the given WCS file; write plot to given PNG filename')
    parser.add_option('--jobid', '-i', dest='solved_id', type=int,help='retrieve result for jobId instead of submitting new image')
    parser.add_option('--substatus', '-s', dest='sub_id', help='Get status of a submission')
    parser.add_option('--jobstatus', '-j', dest='job_id', help='Get status of a job')
    parser.add_option('--jobs', '-J', dest='myjobs', action='store_true', help='Get all my jobs')
    parser.add_option('--jobsbyexacttag', '-T', dest='jobs_by_exact_tag', help='Get a list of jobs associated with a given tag--exact match')
    parser.add_option('--jobsbytag', '-t', dest='jobs_by_tag', help='Get a list of jobs associated with a given tag')
    parser.add_option('--solve-time', '-W', dest='solve_time', default=120., type=float, help='Max time allowed in seconds for image solving.')
    parser.add_option( '--private', '-p',
        dest='public',
        action='store_const',
        const='n',
        default='y',
        help='Hide this submission from other users')
    parser.add_option('--allow_mod_sa','-m',
        dest='allow_mod',
        action='store_const',
        const='sa',
        default='d',
        help='Select license to allow derivative works of submission, but only if shared under same conditions of original license')
    parser.add_option('--no_mod','-M',
        dest='allow_mod',
        action='store_const',
        const='n',
        default='d',
        help='Select license to disallow derivative works of submission')
    parser.add_option('--no_commercial','-c',
        dest='allow_commercial',
        action='store_const',
        const='n',
        default='d',
        help='Select license to disallow commercial use of submission')
    opt,args = parser.parse_args()

    if opt.apikey is None:
        # try the environment
        opt.apikey = os.environ.get('AN_API_KEY', None)
    if opt.apikey is None:
        parser.print_help()
        print()
        print(f'{PREFIX}You must either specify --apikey or set AN_API_KEY')
        sys.exit(-1)

    args = {}
    args['apiurl'] = opt.server
    c = Client(**args)
    c.login(opt.apikey)

    # Print out starting message
    img_name = opt.upload.split("\\")[-1]
    print(f'{PREFIX}(Status = starting) {img_name:11s}')
    
    if opt.upload or opt.upload_url:
        if opt.wcs or opt.kmz or opt.newfits or opt.corr or opt.annotate or opt.calibrate:
            opt.wait = True

        kwargs = dict(
            allow_commercial_use=opt.allow_commercial,
            allow_modifications=opt.allow_mod,
            publicly_visible=opt.public)
        if opt.scale_lower and opt.scale_upper:
            kwargs.update(scale_lower=opt.scale_lower,
                          scale_upper=opt.scale_upper,
                          scale_type='ul')
        elif opt.scale_est and opt.scale_err:
            kwargs.update(scale_est=opt.scale_est,
                          scale_err=opt.scale_err,
                          scale_type='ev')
        elif opt.scale_lower or opt.scale_upper:
            kwargs.update(scale_type='ul')
            if opt.scale_lower:
                kwargs.update(scale_lower=opt.scale_lower)
            if opt.scale_upper:
                kwargs.update(scale_upper=opt.scale_upper)

        for key in ['scale_units', 'center_ra', 'center_dec', 'radius',
                    'downsample_factor', 'positional_error', 'tweak_order', 'crpix_center',
                    'album']:
            if getattr(opt, key) is not None:
                kwargs[key] = getattr(opt, key)
        if opt.parity is not None:
            kwargs.update(parity=int(opt.parity))

        if opt.upload:
            upres = c.upload(opt.upload, **kwargs)

        if opt.upload_url:
            upres = c.url_upload(opt.upload_url, **kwargs)

        stat = upres['status']
        if stat != 'success':
            print(f'{PREFIX}Upload failed: status', stat)
            print(upres)
            sys.exit(-1)

        opt.sub_id = upres['subid']

    if opt.wait:
        if opt.solved_id is None:
            if opt.sub_id is None:
                print(f"{PREFIX}Can't --wait without a submission id or job id!")
                sys.exit(-1)

            while True:
                stat = c.sub_status(opt.sub_id, justdict=True)
                if stat['processing_started'] == 'None':
                    print(f"{PREFIX}(Status = no job  ) {img_name:11s}")

                jobs = stat.get('jobs', [])
                if len(jobs):
                    for j in jobs:
                        if j is not None:
                            break
                    if j is not None:
                        opt.solved_id = j
                        break
                time.sleep(5)

        time_start = time.time()
        time_solving = 0.0
        while True:
            stat = c.job_status(opt.solved_id, justdict=True)
            if stat.get('status','') in ['success']:
                success = (stat['status'] == 'success')
                break
            elif stat.get('status','') in ['failure']:
                print(
                    f"{PREFIX}(Status = {pc.RED}FAILED{pc.END}  ) "
                    f"{img_name:11s} Image solve failed"
                )
                sys.exit(-1)
            elif time_solving > opt.solve_time:
                print(
                    f"{PREFIX}(Status = {pc.RED}FAILED{pc.END}  ) "
                    f"{img_name:11s} Solve time exceeded {opt.solve_time:.0f}s"
                )
                sys.exit(-1)
            else:
                print(f"{PREFIX}(Status = {stat['status']} ) {img_name:11s}")
            time.sleep(5)
            time_solving = time.time() - time_start

    if opt.solved_id:
        # we have a jobId for retrieving results
        retrieveurls = []
        if opt.wcs:
            # We don't need the API for this, just construct URL
            url = opt.server.replace('/api/', '/wcs_file/%i' % opt.solved_id)
            retrieveurls.append((url, opt.wcs))
        if opt.kmz:
            url = opt.server.replace('/api/', '/kml_file/%i/' % opt.solved_id)
            retrieveurls.append((url, opt.kmz))
        if opt.newfits:
            url = opt.server.replace('/api/', '/new_fits_file/%i/' % opt.solved_id)
            retrieveurls.append((url, opt.newfits))
        if opt.corr:
            url = opt.server.replace('/api/', '/corr_file/%i' % opt.solved_id)
            retrieveurls.append((url, opt.corr))

        for url,fn in retrieveurls:
            
            # Try retrieving URL
            delay = 5.0       # Seconds between retry attempts
            max_retries = 50  # Max number of retries
            for attempts in range(max_retries):
                try:
                    f = urlopen(url)
                    txt = f.read()
                    w = open(fn, 'wb')
                    w.write(txt)
                    w.close()
                    break
                except Exception as e:
                    if attempts < max_retries-1:
                        time.sleep(delay)
                    else:
                        print(
                            f'{PREFIX}(Status = {pc.RED}FAILED{pc.END}  ) '
                            f'{img_name:11s} Could not retrieve results'
                        )
                        sys.exit(-1)
            

        if opt.annotate:
            result = c.annotate_data(opt.solved_id)
            with open(opt.annotate,'w') as f:
                f.write(python2json(result))
				
        if opt.calibrate:
            result = c.calibrate_data(opt.solved_id)
            with open(opt.calibrate,'w') as f:
                f.write(python2json(result))

    if opt.wait:
        # behaviour as in old implementation
        opt.sub_id = None

    if opt.sdss_wcs:
        (wcsfn, outfn) = opt.sdss_wcs
        c.sdss_plot(outfn, wcsfn)
    if opt.galex_wcs:
        (wcsfn, outfn) = opt.galex_wcs
        c.galex_plot(outfn, wcsfn)

    if opt.sub_id:
        print(c.sub_status(opt.sub_id))
    if opt.job_id:
        print(c.job_status(opt.job_id))

    if opt.jobs_by_tag:
        tag = opt.jobs_by_tag
        print(c.jobs_by_tag(tag, None))
    if opt.jobs_by_exact_tag:
        tag = opt.jobs_by_exact_tag
        print(c.jobs_by_tag(tag, 'yes'))

    if opt.myjobs:
        jobs = c.myjobs()
        print(jobs)

    print(f"{PREFIX}(Status = {pc.CYAN}finished{pc.END}) {img_name:11s}")