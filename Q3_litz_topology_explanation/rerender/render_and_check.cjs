const { chromium }=require('C:/Users/zq257/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const path=require('path');const fs=require('fs');const {pathToFileURL}=require('url');
(async()=>{
 const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-webgl','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:1120,height:1250},deviceScaleFactor:1.5,colorScheme:'light'});
  const errors=[];page.on('pageerror',e=>errors.push(e.message));
  await page.route('https://**',r=>r.abort()); // verify that the model works offline
  await page.goto(pathToFileURL(path.join(__dirname,'litz-interactive.html')).href,{waitUntil:'networkidle'});
  const f=page.frameLocator('iframe');
  await f.locator('#litz-volume-explainer[data-rendered="bundle"]').waitFor();
  const canvas=f.locator('#litz-solid');
  const before=await canvas.screenshot();
  const box=await canvas.boundingBox();
  await page.mouse.move(box.x+box.width*.5,box.y+box.height*.5);await page.mouse.down();
  await page.mouse.move(box.x+box.width*.65,box.y+box.height*.58,{steps:12});await page.mouse.up();
  await page.waitForTimeout(150);const rotated=await canvas.screenshot();
  if(before.equals(rotated))throw new Error('Mouse rotation did not change the rendered image');
  await f.locator('[data-view="principle"]').click();await f.locator('#litz-next').click();await page.waitForTimeout(550);
  const phase=await f.locator('#litz-phase-value').innerText();if(!phase.includes('1.00'))throw new Error('Step control did not work');
  await f.locator('[data-view="close"]').click();
  await f.locator('#litz-scale').fill('90');await f.locator('#litz-scale').dispatchEvent('input');
  const gap=await f.locator('#litz-gap').innerText();if(!gap.includes('不能保证'))throw new Error('Parameter control did not update the clearance');
  await f.locator('#litz-scale').fill('100');await f.locator('#litz-scale').dispatchEvent('input');
  const frame=page.frames().find(x=>x!==page.mainFrame());
  async function label(title,subtitle,footer){
   await frame.evaluate(({title,subtitle,footer})=>{
    const wrap=document.querySelector('#litz-canvas-wrap');wrap.style.height='760px';wrap.style.background='white';
    wrap.querySelectorAll('.static-label').forEach(x=>x.remove());
    for(const [text,top,bottom,fontSize,weight] of [[title,'10px','auto','25px','500'],[subtitle,'52px','auto','18px','400'],[footer,'auto','12px','19px','400']]){
     const el=document.createElement('div');el.className='static-label';el.textContent=text;
     Object.assign(el.style,{position:'absolute',left:'16px',right:'16px',top,bottom,fontSize,fontWeight:weight,color:'#243541',textAlign:'center',pointerEvents:'none',fontFamily:'Microsoft YaHei, sans-serif'});wrap.append(el);
    }
   },{title,subtitle,footer});await page.waitForTimeout(160);
  }
  await f.locator('[data-view="bundle"]').click();
  await label('484 根有限直径铜管的三维局部','蓝色、橙色追踪两根单丝；其余单丝显示铜表面','轴向局部 4.132 mm  ·  铜径 125.63 μm  ·  各轴等比例');
  await f.locator('#litz-canvas-wrap').screenshot({path:path.join(__dirname,'02-484-copper-tubes.png')});
  await f.locator('[data-view="close"]').click();
  await label('转弯间隙放大：两根单丝及假设漆膜','半透明外层为漆膜；铜和绝缘外径按实际比例显示','假设漆膜 8 μm / 侧  ·  全模型间隙保守下界约 7.54 μm');
  await f.locator('#litz-canvas-wrap').screenshot({path:path.join(__dirname,'03-turning-clearance.png')});
  const result={offline:true,rotation_changed_pixels:true,phase,compression_result:gap,script_errors:errors};
  fs.writeFileSync(path.join(__dirname,'interaction_checks.json'),JSON.stringify(result,null,2));console.log(JSON.stringify(result,null,2));
  if(errors.length)process.exitCode=1;
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
