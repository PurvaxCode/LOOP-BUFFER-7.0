import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

// ─────────────────────────────────────────────
//  MARKER
// ─────────────────────────────────────────────
class Marker {
    String sequence;
    int position;
    String diseaseName;

    Marker(String sequence, int position, String diseaseName) {
        this.sequence    = sequence;
        this.position    = position;
        this.diseaseName = diseaseName;
    }
}

// ─────────────────────────────────────────────
//  DATA LOADER
// ─────────────────────────────────────────────
class DataLoader {
   private static final String HEALTHY_FILE = "healthy_dataset_updated.csv";
private static final String TEST_FILE    = "population members 1 to 300.csv";
private static final String MARKER_FILE  = "disease_markers_updated (1).csv";

    private List<String> healthyPopulation = new ArrayList<>();
    private List<String> testPopulation    = new ArrayList<>();
    private List<Marker>  markers          = new ArrayList<>();

    DataLoader() {
        healthyPopulation = loadDNA(HEALTHY_FILE);
        testPopulation    = loadDNA(TEST_FILE);
        markers           = loadMarkers(MARKER_FILE);
    }

    private List<String> loadDNA(String filePath) {
        List<String> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line; boolean skip = true;
            while ((line = br.readLine()) != null) {
                if (skip) { skip = false; continue; }
                String t = line.trim().toUpperCase();
                if (!t.isEmpty()) list.add(t);
            }
        } catch (Exception e) {
            System.err.println("Could not load: " + filePath + " — " + e.getMessage());
        }
        return list;
    }

    private List<Marker> loadMarkers(String filePath) {
        List<Marker> list = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line; boolean skip = true;
            while ((line = br.readLine()) != null) {
                if (skip) { skip = false; continue; }
                String[] p = line.split(",");
                if (p.length < 3) continue;
                list.add(new Marker(p[1].trim().toUpperCase(),
                                    Integer.parseInt(p[2].trim()),
                                    p[0].trim()));
            }
        } catch (Exception e) {
            System.err.println("Could not load: " + filePath + " — " + e.getMessage());
        }
        return list;
    }

    List<String> getHealthyPopulation() { return healthyPopulation; }
    List<String> getTestPopulation()    { return testPopulation;    }
    List<Marker>  getMarkers()          { return markers;           }
}

// ─────────────────────────────────────────────
//  REFERENCE GENOME BUILDER
// ─────────────────────────────────────────────
class ReferenceGenomeBuilder {
    private static final int K = 4;

    private static int encodeBase(char c) {
        switch (c) { case 'A': return 0; case 'T': return 1; case 'C': return 2; case 'G': return 3; }
        return 0;
    }
    private static int encodeKmer(String dna, int s) {
        int v = 0;
        for (int i = s; i < s + K; i++) { v <<= 2; v |= encodeBase(dna.charAt(i)); }
        return v;
    }
    private static String decodeKmer(int v) {
        char[] m = {'A','T','C','G'}, r = new char[K];
        for (int i = K-1; i >= 0; i--) { r[i] = m[v & 3]; v >>= 2; }
        return new String(r);
    }

    static String buildReference(List<String> pop) {
        if (pop.isEmpty()) return "";
        int len = pop.get(0).length(), pos = len - K + 1;
        Map<Integer, Map<Integer,Integer>> freq = new HashMap<>();
        for (int p = 0; p < pos; p++) freq.put(p, new HashMap<>());
        for (String dna : pop)
            for (int p = 0; p < pos; p++) {
                int enc = encodeKmer(dna, p);
                Map<Integer,Integer> mm = freq.get(p);
                mm.put(enc, mm.getOrDefault(enc, 0) + 1);
            }
        StringBuilder sb = new StringBuilder();
        for (int p = 0; p < pos; p++) {
            int best = 0, max = 0;
            for (Map.Entry<Integer,Integer> e : freq.get(p).entrySet())
                if (e.getValue() > max) { max = e.getValue(); best = e.getKey(); }
            String k = decodeKmer(best);
            if (p == 0) sb.append(k); else sb.append(k.charAt(3));
        }
        return sb.toString();
    }
}

// ─────────────────────────────────────────────
//  DISEASE DETECTOR
// ─────────────────────────────────────────────
class DiseaseDetector {
    private static int mismatches(String a, String b) {
        int c = 0;
        for (int i = 0; i < a.length(); i++) if (a.charAt(i) != b.charAt(i)) c++;
        return c;
    }

    static Map<String,String> detectDiseaseMap(String dna, List<Marker> markers, Set<String> selected) {
        Map<String,String> res = new LinkedHashMap<>();
        for (Marker mk : markers) {
            if (!selected.contains(mk.diseaseName)) continue;
            int pos = mk.position - 1;
            if (pos >= 0 && pos + mk.sequence.length() <= dna.length()) {
                int mm = mismatches(dna.substring(pos, pos + mk.sequence.length()), mk.sequence);
                if      (mm == 0) res.put(mk.diseaseName, "Exact Match");
                else if (mm == 1) res.put(mk.diseaseName, "Mutation");
            }
        }
        return res;
    }
}

// ─────────────────────────────────────────────
//  SIMILARITY CALCULATOR
// ─────────────────────────────────────────────
class SimilarityCalculator {
    static double calculate(String ref, String dna) {
        int n = ref.length(), m = dna.length();
        int[][] dp = new int[n+1][m+1];
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                dp[i][j] = (ref.charAt(i-1) == dna.charAt(j-1))
                    ? dp[i-1][j-1]+1 : Math.max(dp[i-1][j], dp[i][j-1]);
        return (double) dp[n][m] / n * 100.0;
    }
}

// ─────────────────────────────────────────────
//  POINTER CALCULATOR
// ─────────────────────────────────────────────
class PointerCalculator {
    static int similarityPointer(double s) {
        if (s>=85) return 10; if (s>=75) return 9; if (s>=65) return 8;
        if (s>=55) return 7;  if (s>=45) return 6; if (s>=35) return 5;
        if (s>=25) return 4;  if (s>=15) return 3; if (s>=5)  return 2;
        return 1;
    }
    static int diseasePointer(int count) {
        if (count==0) return 10; if (count==1) return 7; return 5;
    }
    static double weightedPointer(int dp, int sp) { return 0.70*dp + 0.30*sp; }

    static String severityLabel(double wp) {
        if (wp >= 9.0) return "Healthy";
        if (wp >= 7.5) return "Monitor Closely";
        if (wp >= 6.0) return "Mild Risk";
        if (wp >= 4.5) return "Moderate Risk";
        if (wp >= 3.0) return "High Risk — Consult Doctor";
        return "Needs Treatment Urgently";
    }
    static Color severityColor(double wp) {
        if (wp >= 9.0) return new Color(0,200,120);
        if (wp >= 7.5) return new Color(60,200,150);
        if (wp >= 6.0) return new Color(220,190,40);
        if (wp >= 4.5) return new Color(240,140,30);
        if (wp >= 3.0) return new Color(230,60,50);
        return new Color(180,0,40);
    }
}

// ─────────────────────────────────────────────
//  MEMBER RESULT
// ─────────────────────────────────────────────
class MemberResult {
    int index;
    double similarity, weightedPointer;
    int simPointer, diseasePointer;
    Map<String,String> diseaseMap;
    String severityLabel;
    Color  severityColor;

    MemberResult(int idx, double sim, Map<String,String> dm) {
        index           = idx;
        similarity      = sim;
        diseaseMap      = dm;
        simPointer      = PointerCalculator.similarityPointer(sim);
        diseasePointer  = PointerCalculator.diseasePointer(dm.size());
        weightedPointer = PointerCalculator.weightedPointer(diseasePointer, simPointer);
        severityLabel   = PointerCalculator.severityLabel(weightedPointer);
        severityColor   = PointerCalculator.severityColor(weightedPointer);
    }
}

// ─────────────────────────────────────────────
//  DNA LOADING PANEL  (helix animation ~5 sec)
// ─────────────────────────────────────────────
class DNALoadingPanel extends JPanel {

    private float  progress = 0f;
    private int    animOff  = 0;
    private javax.swing.Timer swingTimer;

    private static final Color BG     = new Color(8,12,28);
    private static final Color CYAN   = new Color(0,220,255);
    private static final Color GREEN  = new Color(0,255,140);
    private static final Color PURPLE = new Color(160,60,255);
    private static final Color PINK   = new Color(255,80,180);
    private static final Color GOLD   = new Color(255,210,50);
    private static final Color[] RUNG = {GREEN, CYAN, PURPLE, PINK, GOLD};

    DNALoadingPanel(Runnable onDone) {
        setBackground(BG);
        setOpaque(true);
        final int TOTAL = 300;   // ~5 seconds at 16 ms per frame
        final int[] frame = {0};
        swingTimer = new javax.swing.Timer(16, e -> {
            frame[0]++;
            progress = Math.min(1f, (float) frame[0] / TOTAL);
            animOff += 2;
            repaint();
            if (frame[0] >= TOTAL) { swingTimer.stop(); onDone.run(); }
        });
        swingTimer.start();
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D) g.create();
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,      RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        int w = getWidth(), h = getHeight();

        // Background
        g2.setPaint(new GradientPaint(0,0,new Color(5,10,25),0,h,new Color(15,5,40)));
        g2.fillRect(0,0,w,h);

        drawHelix(g2, w, h);

        // Progress bar
        int barW = Math.min(600,(int)(w*0.6)), barH = 18;
        int barX = (w-barW)/2, barY = h-140;

        g2.setColor(new Color(255,255,255,30));
        g2.fillRoundRect(barX, barY, barW, barH, barH, barH);

        int fill = (int)(barW * progress);
        if (fill > 0) {
            g2.setPaint(new GradientPaint(barX,0,CYAN,barX+barW,0,PURPLE));
            g2.fillRoundRect(barX, barY, fill, barH, barH, barH);
            g2.setColor(new Color(0,220,255,50));
            g2.setStroke(new BasicStroke(6));
            g2.drawRoundRect(barX, barY-3, fill, barH+6, barH+6, barH+6);
        }

        // Percentage text
        String pct = (int)(progress*100) + "%";
        g2.setFont(new Font("Courier New", Font.BOLD, 28));
        FontMetrics fm = g2.getFontMetrics();
        g2.setColor(GREEN);
        g2.drawString(pct, (w-fm.stringWidth(pct))/2, barY-22);

        // Title
        g2.setFont(new Font("Courier New", Font.BOLD, 24));
        fm = g2.getFontMetrics();
        g2.setColor(CYAN);
        String title = "SCANNING DNA SEQUENCES...";
        g2.drawString(title, (w-fm.stringWidth(title))/2, barY-66);

        // Subtitle
        g2.setFont(new Font("Courier New", Font.PLAIN, 13));
        fm = g2.getFontMetrics();
        g2.setColor(new Color(160,200,255,180));
        String sub = "Matching markers — please wait";
        g2.drawString(sub, (w-fm.stringWidth(sub))/2, barY+barH+28);

        g2.dispose();
    }

    private void drawHelix(Graphics2D g2, int w, int h) {
        int cx = w/2, hw = 180, sy = 70, ey = h-200, steps = 80;
        int[][] s1 = new int[steps+1][2], s2 = new int[steps+1][2];
        for (int i = 0; i <= steps; i++) {
            double t = (double)i/steps, a = t*4*Math.PI - animOff*0.04;
            int y = sy + (int)(t*(ey-sy));
            s1[i][0] = cx + (int)(Math.cos(a)*hw/2);  s1[i][1] = y;
            s2[i][0] = cx - (int)(Math.cos(a)*hw/2);  s2[i][1] = y;
        }
        // Rungs
        for (int i = 0; i <= steps; i+=5) {
            Color rc = RUNG[(i/5 + animOff/10) % RUNG.length];
            g2.setColor(new Color(rc.getRed(),rc.getGreen(),rc.getBlue(),110));
            g2.setStroke(new BasicStroke(2));
            g2.drawLine(s1[i][0],s1[i][1],s2[i][0],s2[i][1]);
            g2.setColor(rc);
            g2.fillOval(s1[i][0]-5,s1[i][1]-5,10,10);
            g2.fillOval(s2[i][0]-5,s2[i][1]-5,10,10);
        }
        // Strands
        for (int strand = 0; strand < 2; strand++) {
            int[][] pts = (strand==0) ? s1 : s2;
            Color   c   = (strand==0) ? CYAN : GREEN;
            for (int i = 0; i < pts.length-1; i++) {
                int alpha = 80 + (int)(120.0*i/pts.length);
                g2.setColor(new Color(c.getRed(),c.getGreen(),c.getBlue(),alpha));
                g2.setStroke(new BasicStroke(3));
                g2.drawLine(pts[i][0],pts[i][1],pts[i+1][0],pts[i+1][1]);
            }
        }
    }
}

// ─────────────────────────────────────────────
//  DASHBOARD PANEL
// ─────────────────────────────────────────────
class DashboardPanel extends JPanel {

    private static final Color BG      = new Color(8,12,28);
    private static final Color CARD_BG = new Color(18,24,50);
    private static final Color ACCENT  = new Color(0,220,255);
    private static final Color TEXT    = new Color(210,225,255);
    private static final Color MUTED   = new Color(110,130,180);

    private final Set<String> selected = new LinkedHashSet<>();
    private JButton analyzeBtn;

    DashboardPanel(List<String> diseases, Runnable onAnalyze) {
        setBackground(BG);
        setLayout(new BorderLayout());

        // Header
        JPanel header = new JPanel(new BorderLayout());
        header.setBackground(new Color(12,18,40));
        header.setBorder(new EmptyBorder(26,40,20,40));
        JLabel logo = new JLabel("⬡  GENOMICS SCAN");
        logo.setFont(new Font("Courier New", Font.BOLD, 22));
        logo.setForeground(ACCENT);
        JLabel sub = new JLabel("Multi-Disease Genetic Risk Analysis Platform");
        sub.setFont(new Font("Segoe UI", Font.PLAIN, 13));
        sub.setForeground(MUTED);
        JPanel stack = new JPanel(); stack.setOpaque(false);
        stack.setLayout(new BoxLayout(stack, BoxLayout.Y_AXIS));
        stack.add(logo); stack.add(Box.createVerticalStrut(5)); stack.add(sub);
        header.add(stack, BorderLayout.WEST);
        add(header, BorderLayout.NORTH);

        // Card wrapper
        JPanel wrap = new JPanel(new GridBagLayout());
        wrap.setBackground(BG);
        wrap.setBorder(new EmptyBorder(30,60,30,60));

        JPanel card = new JPanel(new BorderLayout(0,20));
        card.setBackground(CARD_BG);
        card.setBorder(new CompoundBorder(
            new LineBorder(new Color(0,220,255,55),1),
            new EmptyBorder(30,36,28,36)));

        JLabel cardTitle = new JLabel("Select Diseases to Analyse");
        cardTitle.setFont(new Font("Segoe UI", Font.BOLD, 19));
        cardTitle.setForeground(TEXT);
        card.add(cardTitle, BorderLayout.NORTH);

        // Checkboxes
        JPanel cbPanel = new JPanel(new GridLayout(0,2,16,12));
        cbPanel.setBackground(CARD_BG);
        for (String d : diseases) {
            JCheckBox cb = new JCheckBox("  " + d);
            cb.setFont(new Font("Segoe UI", Font.PLAIN, 15));
            cb.setForeground(TEXT); cb.setBackground(CARD_BG);
            cb.setFocusPainted(false);
            cb.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
            cb.addItemListener(e -> {
                if (cb.isSelected()) selected.add(d); else selected.remove(d);
                analyzeBtn.setEnabled(!selected.isEmpty());
            });
            cbPanel.add(cb);
        }
        JScrollPane scroll = new JScrollPane(cbPanel);
        scroll.setBorder(BorderFactory.createEmptyBorder());
        scroll.getViewport().setBackground(CARD_BG);
        scroll.setBackground(CARD_BG);
        card.add(scroll, BorderLayout.CENTER);

        // Buttons
        JPanel btnRow = new JPanel(new FlowLayout(FlowLayout.RIGHT,10,0));
        btnRow.setOpaque(false);
        JButton selAll = makeOutlineBtn("Select All");
        selAll.addActionListener(e -> {
            selected.addAll(diseases);
            for (Component c : cbPanel.getComponents())
                if (c instanceof JCheckBox) ((JCheckBox)c).setSelected(true);
            analyzeBtn.setEnabled(true);
        });
        analyzeBtn = makeFilledBtn("  ▶  Run Analysis  ");
        analyzeBtn.setEnabled(false);
        analyzeBtn.addActionListener(e -> onAnalyze.run());
        btnRow.add(selAll); btnRow.add(analyzeBtn);
        card.add(btnRow, BorderLayout.SOUTH);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.fill = GridBagConstraints.BOTH; gbc.weightx = gbc.weighty = 1;
        wrap.add(card, gbc);
        add(wrap, BorderLayout.CENTER);
    }

    Set<String> getSelected() { return selected; }

    private JButton makeFilledBtn(String txt) {
        JButton btn = new JButton(txt) {
            @Override protected void paintComponent(Graphics g) {
                Graphics2D g2 = (Graphics2D) g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                if (isEnabled())
                    g2.setPaint(new GradientPaint(0,0,new Color(0,180,220),getWidth(),0,new Color(100,60,240)));
                else g2.setColor(new Color(50,60,90));
                g2.fillRoundRect(0,0,getWidth(),getHeight(),10,10);
                g2.setColor(isEnabled() ? Color.WHITE : new Color(120,130,160));
                g2.setFont(getFont());
                FontMetrics fm = g2.getFontMetrics();
                g2.drawString(getText(),(getWidth()-fm.stringWidth(getText()))/2,
                    (getHeight()+fm.getAscent()-fm.getDescent())/2);
                g2.dispose();
            }
        };
        btn.setFont(new Font("Segoe UI", Font.BOLD, 14));
        btn.setBorderPainted(false); btn.setContentAreaFilled(false); btn.setFocusPainted(false);
        btn.setPreferredSize(new Dimension(185,42));
        btn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        return btn;
    }

    private JButton makeOutlineBtn(String txt) {
        JButton btn = new JButton(txt) {
            @Override protected void paintComponent(Graphics g) {
                Graphics2D g2 = (Graphics2D) g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setColor(new Color(0,220,255,28)); g2.fillRoundRect(0,0,getWidth(),getHeight(),10,10);
                g2.setColor(new Color(0,220,255,140)); g2.setStroke(new BasicStroke(1.4f));
                g2.drawRoundRect(1,1,getWidth()-2,getHeight()-2,10,10);
                g2.setColor(new Color(0,220,255));
                g2.setFont(getFont()); FontMetrics fm = g2.getFontMetrics();
                g2.drawString(getText(),(getWidth()-fm.stringWidth(getText()))/2,
                    (getHeight()+fm.getAscent()-fm.getDescent())/2);
                g2.dispose();
            }
        };
        btn.setFont(new Font("Segoe UI", Font.PLAIN, 14));
        btn.setBorderPainted(false); btn.setContentAreaFilled(false); btn.setFocusPainted(false);
        btn.setPreferredSize(new Dimension(130,42));
        btn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        return btn;
    }
}

// ─────────────────────────────────────────────
//  RESULTS PANEL
// ─────────────────────────────────────────────
class ResultsPanel extends JPanel {

    private static final Color BG      = new Color(8,12,28);
    private static final Color CARD_BG = new Color(18,24,50);
    private static final Color ACCENT  = new Color(0,220,255);
    private static final Color TEXT    = new Color(210,225,255);
    private static final Color MUTED   = new Color(110,130,180);

    ResultsPanel(List<MemberResult> results, Set<String> selected, Runnable onBack) {
        setBackground(BG);
        setLayout(new BorderLayout(0,0));

        // Header
        JPanel header = new JPanel(new BorderLayout());
        header.setBackground(new Color(12,18,40));
        header.setBorder(new EmptyBorder(18,36,14,36));
        JPanel left = new JPanel(); left.setOpaque(false);
        left.setLayout(new BoxLayout(left,BoxLayout.Y_AXIS));
        JLabel t1 = new JLabel("⬡  ANALYSIS RESULTS");
        t1.setFont(new Font("Courier New",Font.BOLD,20)); t1.setForeground(ACCENT);
        JLabel t2 = new JLabel("Diseases: " + String.join(", ", selected));
        t2.setFont(new Font("Segoe UI",Font.PLAIN,12)); t2.setForeground(MUTED);
        left.add(t1); left.add(Box.createVerticalStrut(4)); left.add(t2);
        header.add(left, BorderLayout.WEST);
        header.add(makeBackBtn(onBack), BorderLayout.EAST);
        add(header, BorderLayout.NORTH);

        // Stats bar + summary cards stacked in CENTER
        JPanel centerStack = new JPanel(new BorderLayout(0, 0));
        centerStack.setBackground(BG);
        centerStack.add(buildStatsBar(results), BorderLayout.NORTH);
        centerStack.add(buildSummaryCards(results), BorderLayout.CENTER);
        add(centerStack,         BorderLayout.CENTER);
        add(buildTable(results), BorderLayout.SOUTH);
    }

    private JPanel buildStatsBar(List<MemberResult> results) {
        int total    = results.size();
        int affected = 0;
        for (MemberResult r : results)
            if (!r.diseaseMap.isEmpty()) affected++;
        double pct = total == 0 ? 0.0 : (affected * 100.0) / total;

        // Two stat tiles in one slim bar
        JPanel bar = new JPanel(new FlowLayout(FlowLayout.CENTER, 30, 0));
        bar.setBackground(BG);
        bar.setBorder(new EmptyBorder(4, 36, 4, 36));

        bar.add(makeStatTile("TOTAL POPULATION", String.valueOf(total),
                new Color(0, 220, 255)));
        bar.add(makeStatTile("AFFECTED MEMBERS", affected + " / " + total,
                new Color(255, 140, 40)));
        bar.add(makeStatTile("AFFECTED %", String.format("%.1f%%", pct),
                pct >= 50 ? new Color(230, 60, 50) : new Color(220, 190, 40)));
        return bar;
    }

    private JPanel makeStatTile(String label, String value, Color accent) {
        JPanel tile = new JPanel(new BorderLayout(0, 4)) {
            @Override protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setColor(new Color(accent.getRed(), accent.getGreen(), accent.getBlue(), 15));
                g2.fillRoundRect(0, 0, getWidth(), getHeight(), 10, 10);
                g2.setColor(new Color(accent.getRed(), accent.getGreen(), accent.getBlue(), 60));
                g2.setStroke(new BasicStroke(1.2f));
                g2.drawRoundRect(0, 0, getWidth()-1, getHeight()-1, 10, 10);
                g2.dispose();
            }
        };
        tile.setOpaque(false);
        tile.setBorder(new EmptyBorder(10, 24, 10, 24));
        tile.setPreferredSize(new Dimension(210, 66));

        JLabel valLbl = new JLabel(value);
        valLbl.setFont(new Font("Courier New", Font.BOLD, 26));
        valLbl.setForeground(accent);
        valLbl.setHorizontalAlignment(SwingConstants.CENTER);

        JLabel lbl = new JLabel(label);
        lbl.setFont(new Font("Segoe UI", Font.PLAIN, 11));
        lbl.setForeground(MUTED);
        lbl.setHorizontalAlignment(SwingConstants.CENTER);

        tile.add(valLbl, BorderLayout.CENTER);
        tile.add(lbl,    BorderLayout.SOUTH);
        return tile;
    }

    private JPanel buildSummaryCards(List<MemberResult> results) {
        int[] buckets = new int[6];
        for (MemberResult r : results) {
            double wp = r.weightedPointer;
            if      (wp >= 9.0) buckets[5]++;
            else if (wp >= 7.5) buckets[4]++;
            else if (wp >= 6.0) buckets[3]++;
            else if (wp >= 4.5) buckets[2]++;
            else if (wp >= 3.0) buckets[1]++;
            else                buckets[0]++;
        }
        String[] labels = {"Urgent Treatment","High Risk","Moderate Risk","Mild Risk","Monitor","Healthy"};
        Color[]  colors = {
            new Color(180,0,40), new Color(230,60,50), new Color(240,140,30),
            new Color(220,190,40), new Color(60,200,150), new Color(0,200,120)
        };
        JPanel grid = new JPanel(new GridLayout(1,6,12,0));
        grid.setBackground(BG);
        grid.setBorder(new EmptyBorder(16,36,12,36));
        for (int i = 5; i >= 0; i--) {
            final int cnt = buckets[i]; final Color c = colors[i];
            JPanel card = new JPanel(new BorderLayout(0,6)) {
                @Override protected void paintComponent(Graphics g) {
                    super.paintComponent(g);
                    Graphics2D g2 = (Graphics2D)g.create();
                    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
                    g2.setColor(new Color(c.getRed(),c.getGreen(),c.getBlue(),18));
                    g2.fillRoundRect(0,0,getWidth(),getHeight(),10,10);
                    g2.setColor(new Color(c.getRed(),c.getGreen(),c.getBlue(),70));
                    g2.setStroke(new BasicStroke(1.2f));
                    g2.drawRoundRect(0,0,getWidth()-1,getHeight()-1,10,10);
                    g2.dispose();
                }
            };
            card.setOpaque(false); card.setBorder(new EmptyBorder(14,16,14,16));
            JLabel num = new JLabel(String.valueOf(cnt));
            num.setFont(new Font("Courier New",Font.BOLD,32)); num.setForeground(c);
            num.setHorizontalAlignment(SwingConstants.CENTER);
            JLabel lbl = new JLabel("<html><center>"+labels[i]+"</center></html>");
            lbl.setFont(new Font("Segoe UI",Font.PLAIN,11)); lbl.setForeground(MUTED);
            lbl.setHorizontalAlignment(SwingConstants.CENTER);
            card.add(num, BorderLayout.CENTER); card.add(lbl, BorderLayout.SOUTH);
            grid.add(card);
        }
        return grid;
    }

    private JPanel buildTable(List<MemberResult> results) {
        String[] cols = {"Member","Diseases Detected","Similarity %","Health Pointer","Severity"};
        Object[][] rows = new Object[results.size()][5];
        for (int i = 0; i < results.size(); i++) {
            MemberResult r = results.get(i);
            StringBuilder sb = new StringBuilder();
            for (Map.Entry<String,String> e : r.diseaseMap.entrySet()) {
                if (sb.length()>0) sb.append("  |  ");
                sb.append(e.getKey()).append(" (").append(e.getValue()).append(")");
            }
            if (sb.length()==0) sb.append("None");
            rows[i][0] = "Member "+r.index;
            rows[i][1] = sb.toString();
            rows[i][2] = String.format("%.1f%%", r.similarity);
            rows[i][3] = String.format("%.2f",   r.weightedPointer);
            rows[i][4] = r.severityLabel;
        }
        DefaultTableModel model = new DefaultTableModel(rows, cols) {
            @Override public boolean isCellEditable(int r, int c) { return false; }
        };
        JTable table = new JTable(model) {
            @Override public Component prepareRenderer(TableCellRenderer rend, int row, int col) {
                Component c = super.prepareRenderer(rend, row, col);
                if (col==4) {
                    c.setBackground(results.get(row).severityColor);
                    c.setForeground(Color.WHITE);
                    if (c instanceof JComponent) ((JComponent)c).setOpaque(true);
                } else {
                    c.setBackground(row%2==0 ? CARD_BG : new Color(22,30,58));
                    c.setForeground(TEXT);
                }
                if (c instanceof JLabel) ((JLabel)c).setBorder(new EmptyBorder(6,12,6,12));
                return c;
            }
        };
        table.setBackground(CARD_BG); table.setForeground(TEXT);
        table.setFont(new Font("Segoe UI",Font.PLAIN,13));
        table.setRowHeight(32);
        table.setGridColor(new Color(28,38,68));
        table.setShowHorizontalLines(true); table.setShowVerticalLines(false);
        table.setIntercellSpacing(new Dimension(0,1));
        table.setSelectionBackground(new Color(0,150,200,80));
        table.setSelectionForeground(TEXT);
        JTableHeader th = table.getTableHeader();
        th.setBackground(new Color(18,26,58)); th.setForeground(ACCENT);
        th.setFont(new Font("Segoe UI",Font.BOLD,13));
        int[] widths = {90,320,90,110,210};
        for (int i=0;i<widths.length;i++) table.getColumnModel().getColumn(i).setPreferredWidth(widths[i]);

        JScrollPane scroll = new JScrollPane(table);
        scroll.setBorder(BorderFactory.createMatteBorder(1,0,0,0,new Color(28,38,68)));
        scroll.getViewport().setBackground(CARD_BG); scroll.setBackground(CARD_BG);

        JPanel wrap = new JPanel(new BorderLayout());
        wrap.setBackground(BG); wrap.setPreferredSize(new Dimension(0,320));
        wrap.add(scroll, BorderLayout.CENTER);
        return wrap;
    }

    private JButton makeBackBtn(Runnable onBack) {
        JButton btn = new JButton("← Back") {
            @Override protected void paintComponent(Graphics g) {
                Graphics2D g2 = (Graphics2D)g.create();
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setColor(new Color(0,220,255,22)); g2.fillRoundRect(0,0,getWidth(),getHeight(),8,8);
                g2.setColor(new Color(0,220,255,110)); g2.setStroke(new BasicStroke(1.2f));
                g2.drawRoundRect(0,0,getWidth()-1,getHeight()-1,8,8);
                g2.setColor(new Color(0,220,255));
                g2.setFont(getFont()); FontMetrics fm = g2.getFontMetrics();
                g2.drawString(getText(),(getWidth()-fm.stringWidth(getText()))/2,
                    (getHeight()+fm.getAscent()-fm.getDescent())/2);
                g2.dispose();
            }
        };
        btn.setFont(new Font("Segoe UI",Font.PLAIN,13));
        btn.setBorderPainted(false); btn.setContentAreaFilled(false); btn.setFocusPainted(false);
        btn.setPreferredSize(new Dimension(100,34));
        btn.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
        btn.addActionListener(e -> onBack.run());
        return btn;
    }
}

// ─────────────────────────────────────────────
//  MAIN  — public class, matches filename exactly
// ─────────────────────────────────────────────
public class GeneticAnalysisApp {

    private JFrame     frame;
    private CardLayout cards;
    private JPanel     root;

    private DataLoader   loader;
    private String       refGenome;
    private List<String> diseaseNames = new ArrayList<>();

    public static void main(String[] args) {
        try { UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName()); }
        catch (Exception ignored) {}
        SwingUtilities.invokeLater(() -> new GeneticAnalysisApp().start());
    }

    private void start() {
        frame = new JFrame("Genomics Scan — Genetic Disease Analysis");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(1100, 740);
        frame.setMinimumSize(new Dimension(900, 600));
        frame.setLocationRelativeTo(null);
        cards = new CardLayout();
        root  = new JPanel(cards);
        frame.setContentPane(root);
        frame.setVisible(true);
        showSplash();
        loadDataAsync();
    }

    private void showSplash() {
        JPanel splash = new JPanel(new GridBagLayout());
        splash.setBackground(new Color(8,12,28));
        JLabel lbl = new JLabel("Loading data...");
        lbl.setFont(new Font("Courier New", Font.BOLD, 20));
        lbl.setForeground(new Color(0,220,255));
        splash.add(lbl);
        root.add(splash, "splash");
        cards.show(root, "splash");
    }

    private void loadDataAsync() {
        new SwingWorker<Void, Void>() {
            @Override protected Void doInBackground() {
                loader    = new DataLoader();
                refGenome = loader.getHealthyPopulation().isEmpty()
                    ? demoRef() : ReferenceGenomeBuilder.buildReference(loader.getHealthyPopulation());
                Set<String> seen = new LinkedHashSet<>();
                for (Marker m : loader.getMarkers()) seen.add(m.diseaseName);
                if (seen.isEmpty()) seen.addAll(Arrays.asList(
                    "Sickle Cell Anemia","Cystic Fibrosis","Huntington's Disease",
                    "BRCA1 Cancer Risk","Tay-Sachs Disease","Phenylketonuria"));
                diseaseNames.addAll(seen);
                return null;
            }
            @Override protected void done() { showDashboard(); }
        }.execute();
    }

    private void showDashboard() {
        DashboardPanel dash = new DashboardPanel(diseaseNames, this::onRunAnalysis);
        root.add(dash, "dashboard");
        cards.show(root, "dashboard");
    }

    private void onRunAnalysis() {
        // Find current dashboard to read selected diseases
        DashboardPanel dash = null;
        for (Component c : root.getComponents())
            if (c instanceof DashboardPanel) { dash = (DashboardPanel) c; break; }
        if (dash == null || dash.getSelected().isEmpty()) return;

        final Set<String> sel = new LinkedHashSet<>(dash.getSelected());

        // Show helix animation; callback triggers analysis
        DNALoadingPanel anim = new DNALoadingPanel(() -> runAnalysisAsync(sel));
        root.add(anim, "loading");
        cards.show(root, "loading");
    }

    private void runAnalysisAsync(Set<String> sel) {
        new SwingWorker<List<MemberResult>, Void>() {
            @Override protected List<MemberResult> doInBackground() {
                List<String> pop     = loader.getTestPopulation();
                List<Marker> markers = loader.getMarkers();
                if (pop.isEmpty()) pop = demoPop(50);
                List<MemberResult> out = new ArrayList<>();
                for (int i = 0; i < pop.size(); i++) {
                    String dna = pop.get(i);
                    double sim = refGenome.isEmpty()
                        ? 55 + Math.random()*35
                        : SimilarityCalculator.calculate(refGenome, dna);
                    Map<String,String> dm = DiseaseDetector.detectDiseaseMap(dna, markers, sel);
                    out.add(new MemberResult(i+1, sim, dm));
                }
                return out;
            }
            @Override protected void done() {
                try {
                    List<MemberResult> res = get();
                    ResultsPanel rp = new ResultsPanel(res, sel, () -> cards.show(root,"dashboard"));
                    root.add(rp, "results");
                    cards.show(root, "results");
                } catch (InterruptedException | ExecutionException ex) { ex.printStackTrace(); }
            }
        }.execute();
    }

    // ── Demo data helpers when CSVs are absent ──
    private String demoRef() {
        StringBuilder sb = new StringBuilder(); String b="ATCG"; Random r=new Random(42);
        for (int i=0;i<500;i++) sb.append(b.charAt(r.nextInt(4))); return sb.toString();
    }
    private List<String> demoPop(int n) {
        List<String> pop=new ArrayList<>(); String b="ATCG"; Random r=new Random();
        for (int i=0;i<n;i++) {
            StringBuilder sb=new StringBuilder();
            for (int j=0;j<500;j++) sb.append(b.charAt(r.nextInt(4))); pop.add(sb.toString());
        }
        return pop;
    }
}
