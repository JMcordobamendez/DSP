# Start time of pyfda_cpp.exe and pyfdax.exe (PyInstaller, PR #4) on Windows: time from the
# start of the process until the main window exists and until it is ready for input
# (its event loop answers a message). The first run after the download is "cold", the following ones "warm".
# Also the memory of the running app (all its processes: pyfdax.exe is a bootloader plus the Python
# process) and the size on disk of both programs, including the folder pyfdax.exe unpacks to %TEMP%.
param([string]$Cpp, [string]$Py, [int]$Runs = 10)

Add-Type @"
using System;
using System.Text;
using System.Collections.Generic;
using System.Runtime.InteropServices;
public static class W32 {
    delegate bool EnumProc(IntPtr h, IntPtr l);
    [DllImport("user32.dll")] static extern bool EnumWindows(EnumProc f, IntPtr l);
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] static extern int GetWindowText(IntPtr h, StringBuilder s, int n);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] public static extern IntPtr SendMessageTimeout(IntPtr h, uint msg, IntPtr w, IntPtr l, uint flags, uint ms, out IntPtr res);
    // visible top level window whose title contains `part`
    public static IntPtr Find(string part) {
        IntPtr found = IntPtr.Zero;
        EnumWindows((h, l) => {
            var sb = new StringBuilder(512);
            GetWindowText(h, sb, 512);
            if (IsWindowVisible(h) && sb.ToString().Contains(part)) { found = h; return false; }
            return true;
        }, IntPtr.Zero);
        return found;
    }
    public static string Titles() {
        var t = new List<string>();
        EnumWindows((h, l) => { var sb = new StringBuilder(512); GetWindowText(h, sb, 512);
                                if (IsWindowVisible(h) && sb.Length > 0) t.Add(sb.ToString()); return true; }, IntPtr.Zero);
        return string.Join(" | ", t);
    }
}
"@

# the main window is found by its title (pyfdax.exe also has a console window)
function Measure-Start([string]$exe, [string]$name, [string]$title, [string[]]$argv) {
    $sw = [Diagnostics.Stopwatch]::StartNew()
    if ($argv) { $p = Start-Process -FilePath $exe -ArgumentList $argv -PassThru } else { $p = Start-Process -FilePath $exe -PassThru }
    $h = [IntPtr]::Zero
    while ($sw.Elapsed.TotalSeconds -lt 60) {
        $h = [W32]::Find($title)
        if ($h -ne [IntPtr]::Zero) { break }
        Start-Sleep -Milliseconds 5
    }
    $t_win = $sw.Elapsed.TotalSeconds
    $t_ready = [double]::NaN
    $found = $h -ne [IntPtr]::Zero
    if ($found) {
        # WM_NULL returns once the event loop of the window processes messages
        $r = [IntPtr]::Zero
        [void][W32]::SendMessageTimeout($h, 0, [IntPtr]::Zero, [IntPtr]::Zero, 0, 60000, [ref]$r)
        $t_ready = $sw.Elapsed.TotalSeconds
    } else {
        Write-Host "not found, visible windows: $([W32]::Titles())"
    }
    # memory once the window has settled, summed over the processes of the app
    Start-Sleep -Milliseconds 2000
    $ps = @(Get-Process -Name $name -ErrorAction SilentlyContinue)
    $ws = ($ps | Measure-Object WorkingSet64 -Sum).Sum / 1MB
    $priv = ($ps | Measure-Object PrivateMemorySize64 -Sum).Sum / 1MB
    $mei = Get-ChildItem $env:TEMP -Directory -Filter '_MEI*' -ErrorAction SilentlyContinue | Sort-Object LastWriteTime | Select-Object -Last 1
    $unpacked = if ($name -eq 'pyfdax' -and $mei) { (Get-ChildItem $mei.FullName -Recurse -File | Measure-Object Length -Sum).Sum / 1MB } else { 0 }
    Get-Process -Name $name -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Milliseconds 1000
    return [pscustomobject]@{ window = $t_win; ready = $t_ready; found = $found; procs = $ps.Count; ws = $ws; priv = $priv; unpacked = $unpacked }
}

$cfg = Join-Path $env:RUNNER_TEMP 'pyfda_cpp_cfg'
$rows = @()
foreach ($app in @(@{ n = 'pyfda_cpp'; exe = $Cpp; t = 'Filter Design Analysis Tool'; a = @('--config-dir', $cfg) },
                   @{ n = 'pyfdax'; exe = $Py; t = 'Filter Design and Analysis'; a = $null })) {
    for ($i = 0; $i -le $Runs; $i++) {
        $r = Measure-Start $app.exe $app.n $app.t $app.a
        Write-Host ("{0} run {1}: window {2:N3} s, ready {3:N3} s, found {4}, {5} processes, working set {6:N1} MB, private {7:N1} MB" -f $app.n, $i, $r.window, $r.ready, $r.found, $r.procs, $r.ws, $r.priv)
        $rows += [pscustomobject]@{ app = $app.n; run = $(if ($i -eq 0) { 'cold' } else { "warm $i" }); window_s = [math]::Round($r.window, 3); ready_s = [math]::Round($r.ready, 3); found = $r.found
                                    procs = $r.procs; working_set_mb = [math]::Round($r.ws, 1); private_mb = [math]::Round($r.priv, 1); unpacked_mb = [math]::Round($r.unpacked, 1) }
    }
}
$rows | Format-Table -AutoSize | Out-String -Width 200 | Write-Output
foreach ($n in 'pyfda_cpp', 'pyfdax') {
    $w = $rows | Where-Object { $_.app -eq $n -and $_.run -ne 'cold' } | Sort-Object ready_s
    $c = $rows | Where-Object { $_.app -eq $n -and $_.run -eq 'cold' }
    $med = $w[[int][math]::Floor($w.Count / 2)].ready_s
    Write-Output ("{0}: cold {1} s, warm median {2} s (min {3}, max {4})" -f $n, $c.ready_s, $med, $w[0].ready_s, $w[-1].ready_s)
    $m = $rows | Where-Object { $_.app -eq $n -and $_.run -ne 'cold' }
    $wsm = ($m | Sort-Object working_set_mb)[[int][math]::Floor($m.Count / 2)].working_set_mb
    $pvm = ($m | Sort-Object private_mb)[[int][math]::Floor($m.Count / 2)].private_mb
    Write-Output ("{0}: memory (warm median) working set {1} MB, private {2} MB" -f $n, $wsm, $pvm)
}

# size on disk
function Size-Of([string]$path) { (Get-ChildItem $path -Recurse -File | Measure-Object Length -Sum).Sum / 1MB }
$cppDir = Split-Path $Cpp
$zip = Join-Path $env:RUNNER_TEMP 'pyfda_cpp.zip'
Compress-Archive -Path "$cppDir\*" -DestinationPath $zip -Force
$un = ($rows | Where-Object { $_.app -eq 'pyfdax' } | Measure-Object unpacked_mb -Maximum).Maximum
Write-Output ("pyfda_cpp: folder {0:N1} MB ({1} files), zipped {2:N1} MB, exe alone {3:N1} MB" -f (Size-Of $cppDir), @(Get-ChildItem $cppDir -Recurse -File).Count, ((Get-Item $zip).Length / 1MB), ((Get-Item $Cpp).Length / 1MB))
Write-Output ("pyfdax: exe {0:N1} MB, unpacked to %TEMP% on every start {1:N1} MB" -f ((Get-Item $Py).Length / 1MB), $un)
$rows | Export-Csv -NoTypeInformation startup_time.csv
