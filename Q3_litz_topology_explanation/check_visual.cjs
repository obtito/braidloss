const { chromium } = require('C:/Users/zq257/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/playwright');
const path = require('path');
(async()=>{
 const browser=await chromium.launch({executablePath:'C:/Program Files/Google/Chrome/Application/chrome.exe',headless:true,args:['--enable-webgl','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 try{
  const page=await browser.newPage({viewport:{width:860,height:1080},deviceScaleFactor:1.5});
  const errors=[];
  page.on('pageerror',e=>errors.push(e.message));
  page.on('console',m=>{if(m.type()==='error')errors.push(m.text());});
  await page.goto('http://127.0.0.1:8874/',{waitUntil:'networkidle'});
  const f=page.frameLocator('iframe');
  await f.locator('#litz-section circle').first().waitFor();
  await f.locator('#litz-volume-explainer').screenshot({path:path.join(__dirname,'01-transposition-principle.png')});
  await f.locator('#litz-next').click();
  await page.waitForTimeout(500);
  const phase=await f.locator('#litz-phase-value').innerText();
  if(!phase.includes('1.00'))throw new Error('Step did not advance: '+phase);
  await f.locator('[data-view="bundle"]').click();
  await f.locator('#litz-volume-explainer[data-rendered="bundle"]').waitFor();
  await f.locator('#litz-volume-explainer').screenshot({path:path.join(__dirname,'02-484-solid-segment.png')});
  await f.locator('[data-view="close"]').click();
  await f.locator('#litz-volume-explainer[data-rendered="close"]').waitFor();
  await f.locator('#litz-volume-explainer').screenshot({path:path.join(__dirname,'03-close-clearance.png')});
  await f.locator('#litz-scale').fill('90');
  await f.locator('#litz-scale').dispatchEvent('input');
  const gap=await f.locator('#litz-gap').innerText();
  if(!gap.includes('不能保证'))throw new Error('Compression warning did not update: '+gap);
  await f.locator('#litz-scale').fill('100');await f.locator('#litz-scale').dispatchEvent('input');
  await page.setViewportSize({width:390,height:1050});
  await f.locator('#litz-volume-explainer').screenshot({path:path.join(__dirname,'04-mobile.png')});
  const frame=page.frames().find(x=>x!==page.mainFrame());
  const bounds=await frame.evaluate(()=>({width:document.documentElement.clientWidth,scroll:document.documentElement.scrollWidth,canvas:document.querySelector('#litz-solid').getBoundingClientRect().width}));
  console.log(JSON.stringify({errors,phase,gap,bounds},null,2));
  if(errors.length)process.exitCode=1;
 } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
