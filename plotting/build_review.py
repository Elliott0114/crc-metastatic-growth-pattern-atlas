"""Print-size caption pages used by the final figure exporters."""
from pathlib import Path
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Bbox
from figure_legends import LEGENDS

OUT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
(OUT/'review/caption_pages').mkdir(parents=True,exist_ok=True)

def wrap(renderer,value,size,width_pt,bold=False):
    font=FontProperties(family='DejaVu Sans',size=size,weight='bold' if bold else 'normal')
    lines=[];current=''
    for word in value.split():
        candidate=(current+' '+word).strip()
        pixels=renderer.get_text_width_height_descent(candidate,font,False)[0]
        if current and pixels>width_pt*100/72:lines.append(current);current=word
        else:current=candidate
    if current:lines.append(current)
    return lines

def caption_page(fig,name):
    """Extend the saved figure below its canvas, retaining its vector artists."""
    title,caption=LEGENDS[name];body_mm=fig._height_mm
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    label=('Supplementary Figure ' if name.startswith('Supplementary') else 'Figure ')+name.rsplit('_',1)[1]+'. '+title
    title_lines=wrap(renderer,label,9,160/25.4*72,True)
    for size in [8.5,8]:
        lines=wrap(renderer,caption,size,160/25.4*72)
        caption_mm=(len(title_lines)*10.5+2+len(lines)*size*1.16)*25.4/72
        if body_mm+caption_mm+3<=225:break
    assert body_mm+caption_mm+3<=225,(name,body_mm,caption_mm)
    total_mm=math.ceil(body_mm+caption_mm+3);added=[];y=1
    for line in title_lines:
        added.append(fig.text(5/170,-y/body_mm,line,fontsize=9,fontweight='bold',va='top'))
        y+=10.5*25.4/72
    y+=2*25.4/72
    for line in lines:
        added.append(fig.text(5/170,-y/body_mm,line,fontsize=size,va='top'))
        y+=size*1.16*25.4/72
    bbox=Bbox.from_bounds(0,-(total_mm-body_mm)/25.4,170/25.4,total_mm/25.4)
    dest=OUT/'review/caption_pages'/name
    fig.savefig(dest.with_suffix('.pdf'),dpi=300,bbox_inches=bbox,pad_inches=0)
    fig.savefig(dest.with_suffix('.png'),dpi=180,bbox_inches=bbox,pad_inches=0)
    for obj in added:obj.remove()
    path=OUT/'review/figure_dimensions.json'
    metrics=json.loads(path.read_text(encoding='utf-8')) if path.exists() else []
    metrics=[m for m in metrics if m['figure']!=name]
    metrics.append({'figure':name,'width_mm':170,'body_height_mm':body_mm,'caption_page_height_mm':total_mm,
                    'caption_font_pt':size,'title_words':len(title.split()),'caption_words':len(caption.split()),
                    'figure_body_export':'original vector artists; only photos/dense spatial layers rasterized'})
    path.write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
