"""Small display-only index of existing GUI SNR metadata; never reads voltages."""
import argparse
import json
from pathlib import Path
from .publish_rvd import atomic_json

def noise_display(model):
    """Project an existing model packet to its displayed day/beam/channels."""
    origin=model['noise_model_origin_unix_s'];stop=model['noise_model_observation_stop_unix_s']
    data=model['noise_model_data'];use=[i for i,t in enumerate(data['times_h']) if origin+t*3600>=stop-86400]
    return dict(schema='pansy_monitor_noise_v1',noise_model_origin_unix_s=origin,
        noise_model_observation_stop_unix_s=stop,noise_zenith_model_status=model.get('noise_zenith_model_status'),
        noise_model_data=dict(times_h=[data['times_h'][i] for i in use],
            noise_power=[[data['noise_power'][i][0]] for i in use]))

def overview(packet):
    first=packet['axis_start_us'];stop=packet['axis_stop_us']
    edges=[];cells=[];column=0
    for frame in sorted(packet['frame_packets'],key=lambda f:f['start_us']):
        left=max(first,frame['start_us']);right=min(stop,frame['stop_us'])
        if right<=left:continue
        if edges and left<edges[-1]:raise ValueError('Overlapping PMSE intervals')
        if edges and left>edges[-1]:edges.append(left);column+=1
        if not edges:edges.append(left)
        peak={}
        for (beam,_,gate),snr in zip(frame['index'],frame['snr_db']):
            if beam==0:peak[gate]=max(peak.get(gate,float('-inf')),snr)
        cells.extend([column,gate,snr] for gate,snr in peak.items())
        edges.append(right);column+=1
    heights=packet['geometry']['altitude_km'][0]
    y=[heights[0]-(heights[1]-heights[0])/2]+[(a+b)/2 for a,b in zip(heights,heights[1:])]+[heights[-1]+(heights[-1]-heights[-2])/2]
    return dict(schema='pansy_zenith_display_v1',selected_day=packet['selected_day'],
        available_days=packet['available_days'],revision=packet['revision'],frame_count=packet['frame_count'],
        latest_mode_us=packet['latest_mode_us'],axis_start_us=first,x=[(t-first)/3.6e9 for t in edges],
        y=y,cells=cells)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('manifest',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--noise-input',type=Path);parser.add_argument('--noise-output',type=Path)
    args=parser.parse_args();atomic_json(args.output,overview(json.loads(args.manifest.read_text())))
    if args.noise_input and args.noise_output:atomic_json(args.noise_output,noise_display(json.loads(args.noise_input.read_text())))

if __name__=='__main__':main()
